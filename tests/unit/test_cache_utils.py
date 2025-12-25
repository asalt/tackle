import json
from pathlib import Path

from tackle.cache import CacheConfig, SnapshotCache


def _write_cache_entry(cache_dir: Path, key_hash: str, last_access: float) -> None:
    data_path = cache_dir / f"cache_{key_hash}.pkl.gz"
    meta_path = cache_dir / f"cache_{key_hash}.json"
    data_path.write_bytes(b"payload")
    meta = {
        "key_hash": key_hash,
        "created_at": last_access,
        "last_access": last_access,
        "schema_version": 1,
    }
    meta_path.write_text(json.dumps(meta))


def test_get_cache_dir_uses_parent_of_results(tmp_path):
    base_dir = tmp_path / "results"

    cache_dir = CacheConfig.default_cache_dir(str(base_dir))

    assert cache_dir == str(tmp_path / ".tackle_cache")


def test_hash_iterable_is_order_independent():
    hash_a, count_a = SnapshotCache.hash_iterable({"b", "a"})
    hash_b, count_b = SnapshotCache.hash_iterable(["a", "b"])

    assert count_a == 2
    assert count_b == 2
    assert hash_a == hash_b


def test_evict_cache_entries_by_count(tmp_path):
    cache_dir = tmp_path
    _write_cache_entry(cache_dir, "a", 1)
    _write_cache_entry(cache_dir, "b", 2)
    _write_cache_entry(cache_dir, "c", 3)

    SnapshotCache.evict_entries(str(cache_dir), max_entries=2, max_bytes=10**9)

    remaining = {p.stem for p in cache_dir.glob("cache_*.json")}
    assert remaining == {"cache_b", "cache_c"}
    assert (cache_dir / "cache_b.pkl.gz").exists()
    assert (cache_dir / "cache_c.pkl.gz").exists()
