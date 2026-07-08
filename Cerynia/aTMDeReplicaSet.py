"""
aTMDeReplicaSet.py unified interface to artemide/Snowflake NP-parameter replica sets.

Merges the old ArtemideReplicaSet and SnowflakeReplicaSet classes from DataProcessor:
both ultimately push parameters into the same harpy.setNPparameters_* functions
(Snowflake's "MAIN" block is exactly artemide's twist-3 "SNOW" block, both landing on
harpy.setNPparameters_tw3), so a single class covers both.

The old positional, version-gated .rep text format is replaced by a self-describing
JSON format:

    {
      "formatVersion": 1,
      "name": "ART25",
      "comment": "...",
      "modules": ["TMDR", "uTMDPDF", "uTMDFF"],
      "initial":  {"TMDR": {"params": [...]}, "uTMDPDF": {"params": [...], "collinearReplica": [0]}, ...},
      "mean":     {...same shape...},
      "replicas": [ {...same shape...}, {...}, ... ]
    }

Each module block carries "params" (a list of floats, or an int replica index for
models that support built-in replica selection) and an optional "collinearReplica"
(a list of integer collinear PDF/FF replica indices, one per hadron).
"""

import json

_NP_SETTERS = {
    "SNOW":         "setNPparameters_tw3",
    "TMDR":         "setNPparameters_TMDR",
    "uTMDPDF":      "setNPparameters_uTMDPDF",
    "uTMDFF":       "setNPparameters_uTMDFF",
    "lpTMDPDF":     "setNPparameters_lpTMDPDF",
    "SiversTMDPDF": "setNPparameters_SiversTMDPDF",
    "wgtTMDPDF":    "setNPparameters_wgtTMDPDF",
}

_COLLINEAR_SETTERS = {
    "uTMDPDF":   "setPDFreplica",
    "uTMDFF":    "setFFreplica",
    "wgtTMDPDF": "sethPDFreplica",
}

MODULES = tuple(_NP_SETTERS)


class aTMDeReplicaSet:

    def __init__(self, name="", comment="", modules=None, initial=None, mean=None, replicas=None):
        modules = list(modules) if modules is not None else []
        unknown = [m for m in modules if m not in _NP_SETTERS]
        if unknown:
            raise ValueError(f"Unknown module(s) {unknown}; must be one of {MODULES}")

        self.name     = name
        self.comment  = comment
        self.modules  = modules
        self.initial  = initial  if initial  is not None else {}
        self.mean     = mean     if mean     is not None else {}
        self.replicas = replicas if replicas is not None else []

    def __repr__(self):
        return (f"<aTMDeReplicaSet: {self.name!r}, {len(self.replicas)} replicas, "
                f"modules={list(self.modules)}>")

    def __len__(self):
        return len(self.replicas)

    @property
    def numberOfReplicas(self):
        return len(self.replicas)

    # --- I/O --------------------------------------------------------------------

    @classmethod
    def from_json(cls, path):
        """Load a replica set from the JSON format described in the module docstring."""
        with open(path) as f:
            data = json.load(f)

        if "modules" not in data:
            raise ValueError(f"'{path}': missing required 'modules' list")

        return cls(
            name     = data.get("name", ""),
            comment  = data.get("comment", ""),
            modules  = data["modules"],
            initial  = data.get("initial", {}),
            mean     = data.get("mean", {}),
            replicas = data.get("replicas", []),
        )

    def to_json(self, path, formatVersion=1):
        """Save the replica set to the JSON format described in the module docstring."""
        data = {
            "formatVersion": formatVersion,
            "name":     self.name,
            "comment":  self.comment,
            "modules":  self.modules,
            "initial":  self.initial,
            "mean":     self.mean,
            "replicas": self.replicas,
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    # --- Access -----------------------------------------------------------------

    def _select(self, num):
        if num == -1:
            return self.initial
        elif num == 0:
            return self.mean
        elif 1 <= num <= len(self.replicas):
            return self.replicas[num - 1]
        else:
            raise IndexError(
                f"Replica number {num} out of range (must be -1, 0, or 1..{len(self.replicas)})"
            )

    def _check_part(self, part):
        if part != "full" and part not in _NP_SETTERS:
            raise ValueError(f"part must be 'full' or one of {MODULES}, got '{part}'")

    def get(self, num=0, part="full"):
        """
        Return the stored data for one replica, without touching artemide.

        num  : -1 = initial replica, 0 = mean replica, 1... = replica from the list
        part : 'full' returns the whole replica dict {module: {"params":..., "collinearReplica":...}};
               a module name (e.g. 'TMDR') returns just that module's block ({} if absent).
        """
        self._check_part(part)
        r = self._select(num)
        return r if part == "full" else r.get(part, {})

    def set(self, num=0, part="full"):
        """
        Send the parameters of one replica to artemide via harpy.

        num  : -1 = initial replica, 0 = mean replica, 1... = replica from the list
        part : 'full' sends every module present in the replica; a module name sends
               only that module.
        """
        import harpy

        self._check_part(part)
        r = self._select(num)
        target_modules = self.modules if part == "full" else [part]

        for module in target_modules:
            block = r.get(module)
            if not block:
                continue

            params = block.get("params")
            if params is not None:
                getattr(harpy, _NP_SETTERS[module])(params)

            collinear = block.get("collinearReplica")
            if collinear:
                setter_name = _COLLINEAR_SETTERS.get(module)
                if setter_name is None:
                    print(f"Modification of {module} collinear replica is not implemented")
                    continue
                setter = getattr(harpy, setter_name)
                for h, n in enumerate(collinear, start=1):
                    setter(int(n), h=h)
