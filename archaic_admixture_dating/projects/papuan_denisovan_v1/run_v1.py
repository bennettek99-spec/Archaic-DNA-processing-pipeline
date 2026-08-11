"""Convenience launcher for the Papuan Denisovan V1 workflow."""

from __future__ import annotations

from archaic_admixture_dating.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["run-all", *(__import__("sys").argv[1:])]))
