import importlib
from CutFinder.configs import ConfigObj, ConfigRef, GlobalConf


class ConfigReader:
    def __init__(self, path):
        filename = path.split("/")[-1].replace(".py", "")
        spec = importlib.util.spec_from_file_location(filename, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        self.glob = None
        self.refs = []
        self.objs = []

        for attr in dir(module):
            obj = getattr(module, attr)
            if isinstance(obj, ConfigObj):
                obj.name = attr
                self.objs.append(obj)
            elif isinstance(obj, ConfigRef):
                obj.name = attr
                self.refs.append(obj)
            elif isinstance(obj, GlobalConf):
                if self.glob is not None:
                    raise ValueError(
                        "Multiple GlobalConf instances found in the config file."
                    )
                self.glob = obj
        if self.glob is None:
            raise ValueError("No GlobalConf instance found in the config file.")
