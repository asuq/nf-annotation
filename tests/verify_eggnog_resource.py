#!/usr/bin/env python3
"""Validate preparation against the actual pinned eggNOG API inside its image."""

import array
import pickle
import sqlite3
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_resources import AnnotationResourceError, prepare_eggnog, run
from eggnogmapper.annotator.e7.annotate import AnnotationEngine
from eggnogmapper.annotator.e7.db import EggnogDB


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "source"
        source.mkdir()
        with sqlite3.connect(source / "eggnog.db") as db:
            db.executescript("""
                CREATE TABLE version(version TEXT);
                INSERT INTO version VALUES ('7.0.0');
                CREATE TABLE protein_names(id INTEGER PRIMARY KEY, name TEXT, taxid INTEGER);
                INSERT INTO protein_names VALUES (1, '2.control', 2);
                CREATE TABLE prots(id INTEGER PRIMARY KEY, gos TEXT, kegg_ko TEXT);
                INSERT INTO prots VALUES (1, 'GO:0003674', 'K00001');
                CREATE TABLE event_index(protein_id INTEGER, events TEXT);
                INSERT INTO event_index VALUES (1, '[]');
                CREATE TABLE sp_events(i INTEGER, name TEXT, og TEXT, og_lca TEXT,
                    ev_lca TEXT, sp_overlap REAL, side1 TEXT, side2 TEXT);
                INSERT INTO sp_events VALUES (1, 'control', 'COG0001', '2', '2', 1, '[1]', '[1]');
                CREATE TABLE ogs(og TEXT);
                INSERT INTO ogs VALUES ('COG0001');
            """)
        with sqlite3.connect(source / "eggnog.taxa.db") as db:
            db.execute("CREATE TABLE fixture(id INTEGER)")
        (source / "eggnog.taxa.db.traverse.pkl").write_bytes(pickle.dumps({}))
        obo = "format-version: 1.2\n"
        for go, namespace in (
            ("0003674", "molecular_function"),
            ("0008150", "biological_process"),
            ("0005575", "cellular_component"),
        ):
            obo += f"\n[Term]\nid: GO:{go}\nname: {namespace}\nnamespace: {namespace}\n"
        (source / "go-basic.obo").write_text(obo)
        with EggnogDB(str(source / "eggnog.db"), load_taxids=False) as db:
            assert db.get_version() == "7.0.0"
            engine = AnnotationEngine(db, go_obo_path=str(source / "go-basic.obo"))
            masks = engine.build_field_presence(
                [dict(row) for row in db.conn.execute("SELECT * FROM prots")], 2
            )
        (source / "eggnog.db.fieldpresence.bin").write_bytes(masks.tobytes())
        (source / "eggnog.db.taxids.bin").write_bytes(
            array.array("i", [0, 2]).tobytes()
        )
        (source / "control.faa").write_text(">2.control\nMALWMRLLPLLALLALWGPDPAAA\n")
        run(
            [
                "diamond",
                "makedb",
                "--in",
                str(source / "control.faa"),
                "--db",
                str(source / "eggnog_proteins.dmnd"),
            ],
            source,
            "makedb",
        )
        destination = root / "prepared"
        destination.mkdir()
        result = prepare_eggnog(source, destination)
        assert result["database_version"] == "7.0.0"
        # A same-size but incorrect cache must fail even with a valid SQLite DB.
        (source / "eggnog.db.taxids.bin").write_bytes(
            array.array("i", [0, 3]).tobytes()
        )
        invalid = root / "invalid"
        invalid.mkdir()
        try:
            prepare_eggnog(source, invalid)
        except AnnotationResourceError as error:
            assert "Taxid cache disagrees" in str(error)
        else:
            raise AssertionError("A mismatched taxid cache was accepted")
    print("Pinned eggNOG API and complete cache comparison controls passed")


if __name__ == "__main__":
    main()
