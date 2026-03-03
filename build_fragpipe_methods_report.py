#!/usr/bin/env python3
"""Compatibility wrapper for the formalized FragPipe search report module."""

from __future__ import annotations

import sys

from tackle.search_report import main


if __name__ == "__main__":
    sys.exit(main())
