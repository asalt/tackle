import pandas as pd

from tackle.cache import CacheConfig, DataSnapshot, DataSnapshotCache


class FakeData:
    def __init__(self, *, experiment_file: str, data_dir: str):
        self.experiment_file = experiment_file
        self.data_dir = data_dir

        self.additional_info = None
        self.geneid_subset = {"1", "2"}
        self.ignore_geneid_subset = None
        self.genefile_norm = None
        self.funcats = None
        self.funcats_inverse = None
        self.cluster_annotate_cols = None
        self.taxon = "all"
        self.unique_pepts = 0
        self.non_zeros = 0
        self.nonzero_subgroup = None
        self.SRA = "S"
        self.number_sra = 1
        self.normalize_across_species = False
        self.export_all = False
        self.metrics = False
        self.metrics_after_filter = True
        self.metrics_unnormed_area = True
        self.ifot = False
        self.ifot_ki = False
        self.ifot_tf = False
        self.median = False
        self.quantile75 = False
        self.quantile90 = False
        self.impute_missing_values = False
        self.fill_na_zero = True
        self.imputation_backend = "gaussian"
        self.lupine_mode = "local"
        self.batch = None
        self.batch_nonparametric = False
        self.batch_noimputation = False
        self.covariate = None
        self.norm_info = None

        self.data = pd.DataFrame({"a": [1]})
        self.df_filtered = pd.DataFrame({"b": [2]})
        self.col_metadata = pd.DataFrame({"group": ["x"]}, index=["sample"])
        self.taxon_ratios = pd.DataFrame({"9606": [0.5]}, index=["sample"])
        self.gid_funcat_mapping = {"1": "A"}
        self._metric_values = {"exp": {"PSMs": {"Total": 1}}}
        self._areas = pd.DataFrame({"sample": [1.0]}, index=["1"])
        self._areas_log = pd.DataFrame({"sample": [0.0]}, index=["1"])
        self._mask = pd.DataFrame({"sample": [False]}, index=["1"])
        self._zeros = pd.DataFrame({"sample": [False]}, index=["1"])
        self.minval = 1.0
        self.batch_applied = None
        self.normed = False

        self._gid_symbol = "present"
        self.panel = "panel"
        self.exps = "exps"


def test_data_snapshot_roundtrip_resets_ephemeral(tmp_path):
    config_file = tmp_path / "config.ini"
    config_file.write_text("[x]\n")
    fake = FakeData(experiment_file=str(config_file), data_dir=str(tmp_path / "data"))

    payload = DataSnapshot.from_data(fake).to_payload()
    restored = FakeData(experiment_file=str(config_file), data_dir=str(tmp_path / "data"))
    restored.data = pd.DataFrame({"a": [999]})

    DataSnapshot.from_payload(payload).apply_to(restored)

    assert restored.data.equals(fake.data)
    assert restored._gid_symbol is None
    assert restored.panel is None
    assert restored.exps is None


def test_data_snapshot_cache_save_and_restore(tmp_path):
    config_file = tmp_path / "config.ini"
    config_file.write_text("[x]\n")

    config = CacheConfig(
        enabled=True,
        cache_dir=str(tmp_path),
        max_entries=5,
        max_bytes=10**9,
        bust=None,
    )

    fake1 = FakeData(experiment_file=str(config_file), data_dir=str(tmp_path / "data"))
    cache1 = DataSnapshotCache(config, only_local=False)
    assert cache1.try_restore(fake1) is False
    assert cache1.maybe_save(fake1) is True

    fake2 = FakeData(experiment_file=str(config_file), data_dir=str(tmp_path / "data"))
    fake2.data = pd.DataFrame({"a": [999]})
    cache2 = DataSnapshotCache(config, only_local=False)
    assert cache2.try_restore(fake2) is True
    assert fake2.data.equals(fake1.data)

    cache3 = DataSnapshotCache(config, only_local=True)
    assert cache3.try_restore(
        FakeData(experiment_file=str(config_file), data_dir=str(tmp_path / "data"))
    ) is False
