from __future__ import annotations

import pytest

from archaic_admixture_dating.checkpointing import CheckpointStore, Deadline, atomic_write_text


def test_completed_unit_is_skipped_only_while_output_is_valid(tmp_path):
    output = tmp_path / "result.txt"
    atomic_write_text(output, "valid\n")
    store = CheckpointStore(tmp_path / "checkpoint.json", "fixture", "abc", 7)
    state = store.load()
    store.mark_completed(state, "chr21", [output])
    assert store.unit_valid(state, "chr21")
    atomic_write_text(output, "changed\n")
    assert not store.unit_valid(state, "chr21")


def test_checkpoint_refuses_changed_configuration(tmp_path):
    first = CheckpointStore(tmp_path / "checkpoint.json", "fixture", "abc", 7)
    first.save(first.load())
    second = CheckpointStore(tmp_path / "checkpoint.json", "fixture", "different", 7)
    with pytest.raises(ValueError, match="configuration mismatch"):
        second.load()


def test_deadline_validates_stop_buffer():
    with pytest.raises(ValueError):
        Deadline(5, 5)
    assert not Deadline(10, 1).should_stop()
