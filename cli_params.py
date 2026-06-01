from __future__ import annotations
from argparse import ArgumentParser
from dataclasses import asdict

class CLIParams:
    def cli(self) -> CLIParams:
        parser = ArgumentParser()
        for k, v in asdict(self).items():
            parser.add_argument(f"--{k}", type=type(v), default=v)
        args = parser.parse_args()

        for k, v in vars(args).items():
            setattr(self, k, v)
        return self

