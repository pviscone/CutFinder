"""
Multithread-safe extraction of RDataFrame columns into flat numpy float32 arrays.

ForeachSlot fills per-slot C++ buffers, so worker threads never share mutable
state and never re-enter Python. With ROOT.EnableImplicitMT() events are
processed in arbitrary order, but all columns extracted in the SAME call are
shuffled identically: every event lands at the same index in every column.
Separate calls are NOT guaranteed to share the same event order.

Scalar columns -> one 1-D array of length n_entries per column.
RVec columns   -> one 1-D array per column with the elements of all entries
                  concatenated; event boundaries are lost, but event chunks
                  appear in the same order in every column of the same call.
"""
import itertools
from contextlib import contextmanager

import numpy as np
import ROOT


# Interpreter-owned scratch buffers keyed by a unique name per extraction:
# one inner vector per slot (scalar path) or per (slot, column) (RVec path).
# _filled_buffers erases the entry after readback, so repeated extractions
# do not accumulate in g_bufvecs.
_CPP_HELPERS = r"""
#include <map>
#include <string>
#include <vector>
namespace Extractor {
std::map<std::string, std::vector<std::vector<float>>> g_bufvecs;
std::vector<std::vector<float>>& get_bufvec(const std::string& name) { return g_bufvecs[name]; }
void init_bufvec(const std::string& name, unsigned int nbufs) { g_bufvecs[name].assign(nbufs, {}); }
void free_bufvec(const std::string& name) { g_bufvecs.erase(name); }
}
"""

_helpers_declared = False


def _declare_helpers():
    """Declare the C++ buffer helpers once per process."""
    global _helpers_declared
    if not _helpers_declared:
        ROOT.gInterpreter.Declare(_CPP_HELPERS)
        _helpers_declared = True


_functor_ids = itertools.count()


def _declare_functor(types, push_lines):
    """JIT-compile a ForeachSlot functor; return (functor, unique C++ name).

    types: GetColumnType() per column. push_lines: one C++ statement per
    column, where `slot` is the worker slot, f<i> the i-th column value and
    `bufvec` the pointer to the interpreter-owned buffers.

    The class name must be unique: cling cannot redefine an already declared
    class, and concurrent extractions must not clobber each other. ForeachSlot
    copies the functor per worker, but every copy writes through the same
    bufvec pointer into disjoint slot buffers, so no mutable state is shared.
    """
    _declare_helpers()
    name = f"extractor_{next(_functor_ids)}"
    params = ", ".join(f"const {t}& f{i}" for i, t in enumerate(types))
    body = "\n".join(" " * 16 + line for line in push_lines)
    ROOT.gInterpreter.Declare(f"""
        namespace Extractor {{
        class {name} {{
        public:
            std::vector<std::vector<float>>* bufvec = nullptr;
            void set_bufvec(std::vector<std::vector<float>>* bv) {{ bufvec = bv; }}
            void operator() (unsigned int slot, {params}) {{
{body}
            }}
        }};
        }}
        """)
    return getattr(ROOT.Extractor, name)(), name


@contextmanager
def _filled_buffers(rdf, columns, types, nbufs, push_lines):
    """Run ForeachSlot with a fresh functor; yield the filled buffers.

    The buffers are erased from the interpreter on exit, so everything needed
    must be copied into numpy inside the with-block.
    """
    extractor, name = _declare_functor(types, push_lines)
    bufname = f"buf_{name}"
    ROOT.Extractor.init_bufvec(bufname, nbufs)
    bufvec = ROOT.Extractor.get_bufvec(bufname)
    extractor.set_bufvec(bufvec)
    try:
        rdf.ForeachSlot(extractor, columns)
        yield bufvec
    finally:
        ROOT.Extractor.free_bufvec(bufname)


def extract_columns(rdf, columns):
    """Extract columns to a dict {name: flat numpy float32 array}.

    `columns` must be all scalar or all RVec; call once per group to mix.
    Multithread requires ROOT.EnableImplicitMT(); works single-threaded too.
    """
    if not columns:
        return {}
    types = [rdf.GetColumnType(c) for c in columns]
    is_rvec = ["RVec" in t for t in types]
    if any(is_rvec) and not all(is_rvec):
        raise TypeError(
            "extract_columns cannot mix scalar and RVec columns; call it once per group"
        )
    nslots = max(ROOT.GetThreadPoolSize(), 1)
    extract = _extract_rvec if all(is_rvec) else _extract_scalar
    return extract(rdf, columns, types, nslots)


def _extract_scalar(rdf, columns, types, nslots):
    """Scalar columns: one 1-D float32 array of length n_entries per column."""
    n_columns = len(columns)
    # Each event's values are pushed contiguously: [e0_c0, e0_c1, ..., e1_c0, ...].
    # Reshaping to (n, n_columns) and slicing mat[:, j] recovers column j, so the
    # arbitrary multithreaded event order stays identical across columns.
    pushes = [
        f"(*bufvec)[slot].push_back(static_cast<float>(f{i}));"
        for i in range(n_columns)
    ]
    with _filled_buffers(rdf, columns, types, nslots, pushes) as bufvec:
        slot_counts = [bufvec[s].size() // n_columns for s in range(nslots)]
        out = {c: np.empty(sum(slot_counts), dtype=np.float32) for c in columns}
        offset = 0
        for s, n in enumerate(slot_counts):
            if n:
                mat = np.asarray(bufvec[s], dtype=np.float32).reshape(n, n_columns)
                for j, c in enumerate(columns):
                    out[c][offset : offset + n] = mat[:, j]
                offset += n
        return out


def _extract_rvec(rdf, columns, types, nslots):
    """RVec columns: one flat 1-D float32 array per column (events concatenated)."""
    n_columns = len(columns)
    # One flat buffer per (slot, column), indexed slot * n_columns + col. Within
    # a slot, events are appended to every column's buffer in the same order,
    # so concatenating slots in order keeps event chunks aligned across columns.
    # n_columns is baked into the generated code, so the functor needs no state.
    pushes = [
        f"for (auto v : f{i}) "
        f"(*bufvec)[slot * {n_columns} + {i}].push_back(static_cast<float>(v));"
        for i in range(n_columns)
    ]
    with _filled_buffers(rdf, columns, types, nslots * n_columns, pushes) as bufvec:
        chunks = {c: [] for c in columns}
        for s in range(nslots):
            for j, c in enumerate(columns):
                v = bufvec[s * n_columns + j]
                if v.size():
                    chunks[c].append(np.asarray(v, dtype=np.float32))
        return {
            c: np.concatenate(chunks[c]) if chunks[c] else np.empty(0, dtype=np.float32)
            for c in columns
        }
