"""Tests for the gauge utilities in clawpack.pyclaw.gauges"""

import numpy as np
import pytest

from clawpack.pyclaw import gauges


def _write_gauge(path, gauge_id, num_eqn=2):
    gauge = gauges.GaugeSolution()
    gauge.id = gauge_id
    gauge.location = (0.0, 0.0)
    gauge.level = np.ones(3, dtype=int)
    gauge.t = np.linspace(0.0, 1.0, 3)
    gauge.q = np.arange(num_eqn * 3, dtype=float).reshape((num_eqn, 3))
    gauge.write(str(path))


def test_find_gauge_ids(tmp_path):
    # 5-digit padded ids, an id longer than 5 digits, and non-gauge files
    for gauge_id in (12, 1, 1234567):
        _write_gauge(tmp_path, gauge_id)
    (tmp_path / "gauge_notes.txt").write_text("not a gauge\n")
    (tmp_path / "gauges.data").write_text("not a gauge\n")

    assert gauges.find_gauge_ids(tmp_path) == [1, 12, 1234567]

    # Long ids must round trip through the reader with the found id
    gauge = gauges.GaugeSolution(gauge_id=1234567, path=str(tmp_path))
    assert gauge.id == 1234567


def test_compare_gauges_bad_fields(tmp_path):
    pytest.importorskip("matplotlib")
    for name, num_eqn in (("a", 2), ("b", 3)):
        (tmp_path / name).mkdir()
        _write_gauge(tmp_path / name, 1, num_eqn=num_eqn)
    paths = [str(tmp_path / "a"), str(tmp_path / "b")]

    # Field 2 exists only in the second path
    with pytest.raises(ValueError, match="out of range"):
        gauges.compare_gauges(paths, 1, fields=[0, 2])
