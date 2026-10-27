# Rebuild the v0.1.0 release payloads

Run from the repository root:

```powershell
python scripts/build_public_release_v1.py
```

The default output is the ignored directory
`docs/public_release_v1/local_private/artifacts_v0.1.0/`. The directory must not
already exist. For an independent repeat, choose a new output directory:

```powershell
python scripts/build_public_release_v1.py --out docs/public_release_v1/local_private/artifacts_v0.1.0_repeat
```

The builder uses only the prepared local export, these package notices and the
Python standard library. It does not execute training, prediction, scoring,
random splitting, network downloads or serialized estimators. It preserves all
source files. Archive entries are sorted and use a fixed timestamp, so
independent builds of unchanged inputs produce identical ZIP bytes.

Outputs are `reo-model-assets-v0.1.0.zip`, `reo-metadata-assets-v0.1.0.zip`,
`SHA256SUMS`, `release_file_manifest.json`, `projection_membership_audit.json`,
`sensitive_pattern_audit.json`, `model_hash_preservation.json`,
`public_release_summary.json` and `asset_catalog.json`. The summary and catalog
are safe publication metadata: they describe files, hashes, counts and column
names, without identity rows, fitted vectors or observations. Reviewed copies
may be committed under this directory; ZIPs are release attachments.

Every read is selected from a fixed allowlist. All model bytes are checked
against both the prepared asset inventory and model export manifest. Metadata
projection counts and retained strings are checked, forbidden source fields
are excluded, and four fold hashes are preserved. Source file sizes and
modification times are compared before and after the build. The output path
must remain inside this directory's `local_private` tree, and existing outputs
are never overwritten. Package notices describe the author-owned-rights grant
and upstream materials that are not relicensed.
