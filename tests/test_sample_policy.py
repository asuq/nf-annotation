"""Explicit domain routing and recorded BUSCO-selection contracts."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import sample_policy as policy


class SamplePolicyTests(unittest.TestCase):
    def config(self):
        return {
            "domain": "metadata",
            "busco_domain_lineages": {
                "Bacteria": ["bacteria_odb12", "bacillota_odb12"],
                "Archaea": ["archaea_odb12"],
            },
        }

    def test_explicit_domains_route_to_their_own_primary_lineage(self):
        for name, kingdom in policy.DOMAINS.items():
            row = policy.resolve_sample({"Domain": name}, self.config())
            self.assertEqual(row["annotation_domain"], name)
            self.assertEqual(row["prokka_kingdom"], name)
            self.assertEqual(row["barrnap_kingdom"], kingdom)
            self.assertEqual(
                policy.primary_column(row, "BUSCO_wrong"),
                "BUSCO_" + self.config()["busco_domain_lineages"][name][0],
            )

    def test_legacy_bacterial_policy_is_unchanged(self):
        self.assertEqual(policy.resolve_sample({}, {"domain": "bacteria"}), {})
        self.assertEqual(
            policy.primary_column({}, "BUSCO_bacillota_odb12"), "BUSCO_bacillota_odb12"
        )

    def test_missing_unknown_or_conflicting_domains_fail(self):
        for value in (None, "", "NA", "unknown", "Eukaryota", "Bacteria;Archaea"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                policy.resolve_sample({"Domain": value}, self.config())
        for changed in (
            {"prokka_kingdom": "Bacteria"},
            {"barrnap_kingdom": "bac"},
            {"busco_domain_lineages": {}},
            {"busco_primary_column": "BUSCO_bacteria_odb12"},
        ):
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                policy.resolve_sample({"Domain": "Archaea"}, self.config() | changed)

    def test_profiles_and_published_records_fail_closed(self):
        for value in (
            {"Archaea": ["bacteria_odb12"]},
            {"Bacteria": ["bacteria_odb12"], "Archaea": ["bacteria_odb12"]},
            {"Archaea": []},
            {"Archaea": ["../archaea"]},
        ):
            with self.subTest(value=value), self.assertRaises(ValueError):
                policy.profiles(value)
        row = policy.resolve_sample({"Domain": "Archaea"}, self.config())
        for changed in (
            {"primary_busco_lineage": "bacteria_odb12"},
            {"prokka_kingdom": "Bacteria"},
        ):
            with self.assertRaises(ValueError):
                policy.primary_column(row | changed, None)
        del row["barrnap_kingdom"]
        with self.assertRaises(ValueError):
            policy.primary_column(row, None)
        with self.assertRaises(ValueError):
            policy.qc_scope({"qc_scope": "unknown"})

    def test_domain_is_required_even_with_profiles(self):
        for value in (None, "", "auto", "Bacteria", True):
            with (
                self.subTest(value=value),
                self.assertRaisesRegex(ValueError, "--domain is required"),
            ):
                policy.resolve_sample({}, self.config() | {"domain": value})
        with self.assertRaisesRegex(ValueError, "--domain is required"):
            policy.resolve_sample({}, {})

    def test_single_domain_modes_validate_metadata_and_record_source(self):
        for name in policy.DOMAINS:
            config = self.config() | {"domain": name.lower()}
            for metadata in ({}, {"Domain": ""}, {"Domain": name.upper()}):
                row = policy.resolve_sample(metadata, config)
                self.assertEqual(row["annotation_domain"], name)
                self.assertEqual(row["domain_source"], "domain")
            other = "Archaea" if name == "Bacteria" else "Bacteria"
            with self.assertRaisesRegex(ValueError, "conflicts with --domain"):
                policy.resolve_sample({"Domain": other}, config)
        with self.assertRaisesRegex(ValueError, "conflicts with --domain"):
            policy.resolve_sample({"Domain": "Archaea"}, {"domain": "bacteria"})

    def test_metadata_requires_exact_header_and_explicit_profiles(self):
        with self.assertRaises(ValueError):
            policy.resolve_sample({"domain": "Bacteria"}, self.config())
        for mode in ("metadata", "archaea"):
            with self.assertRaisesRegex(ValueError, "Missing explicit BUSCO"):
                policy.resolve_sample({"Domain": "Archaea"}, {"domain": mode})
        row = policy.resolve_sample({"Domain": "archaea"}, self.config())
        self.assertEqual(row["domain_source"], "metadata:Domain")

    def test_removed_options_and_withdrawn_source_are_not_aliases(self):
        for name in ("domain_column", "prokka_kingdom", "barrnap_kingdom"):
            with (
                self.subTest(name=name),
                self.assertRaisesRegex(ValueError, "Unsupported options"),
            ):
                policy.resolve_sample({}, self.config() | {name: None})
        row = policy.resolve_sample({"Domain": "Archaea"}, self.config())
        with self.assertRaisesRegex(ValueError, "unsupported domain declaration"):
            policy.primary_column(row | {"domain_source": "prokka_kingdom"}, None)
