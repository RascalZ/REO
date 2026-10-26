# Delete only the three author-approved obsolete reading folders after backup.
$ErrorActionPreference = 'Stop'
$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$auditDirectory = Join-Path $repoRoot 'docs/repository_release_v1'
$plan = Get-Content -LiteralPath (Join-Path $auditDirectory 'cleanup_plan.json') -Raw | ConvertFrom-Json
$package = Get-Content -LiteralPath (Join-Path $auditDirectory 'submission_package_manifest.json') -Raw | ConvertFrom-Json
$archivePath = [IO.Path]::GetFullPath((Join-Path $repoRoot $plan.archive))
$packagePath = [IO.Path]::GetFullPath((Join-Path $repoRoot $package.package))
if (-not $plan.archive_CRC_and_hashes_verified -or -not $plan.v4_supplements_migrated_and_verified) { throw 'Backup or migration not verified' }
if ((Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash -ne $plan.archive_sha256) { throw 'Backup digest differs' }
if ((Get-FileHash -LiteralPath $packagePath -Algorithm SHA256).Hash -ne $package.sha256) { throw 'Submission package digest differs' }
$allowedNames = @('author_review_v1', 'author_review_v2', 'author_review_v2_zh')
$checkedTargets = @()
foreach ($target in $plan.targets) {
    $absolute = [IO.Path]::GetFullPath((Join-Path $repoRoot $target.directory))
    $leaf = Split-Path -Leaf $absolute
    $expected = [IO.Path]::GetFullPath((Join-Path (Join-Path $repoRoot 'manuscript_v1') $leaf))
    if ($leaf -notin $allowedNames -or -not $absolute.Equals($expected, [StringComparison]::OrdinalIgnoreCase)) { throw 'Unexpected cleanup target' }
    if (Test-Path -LiteralPath $absolute) {
        $item = Get-Item -LiteralPath $absolute
        if (-not $item.PSIsContainer -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) { throw 'Unexpected target type or link' }
        $children = @(Get-ChildItem -LiteralPath $absolute -Force -Recurse)
        if (@($children | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }).Count) { throw 'Refuse to traverse a reparse point' }
        if (@($children | Where-Object { $_.Name.StartsWith('~$') }).Count) { throw 'An old document is open' }
        foreach ($child in @($children | Where-Object { -not $_.PSIsContainer })) {
            $relative = [IO.Path]::GetRelativePath($repoRoot, $child.FullName).Replace('\', '/')
            if (-not @($plan.all_files_before_cleanup | Where-Object { $_.path -eq $relative }).Count) { throw 'An unrecorded file appeared in an old folder' }
        }
    }
    $checkedTargets += $absolute
}
if ($checkedTargets.Count -ne 3) { throw 'Expected exactly three known targets' }
foreach ($file in $plan.all_files_before_cleanup) {
    $absolute = [IO.Path]::GetFullPath((Join-Path $repoRoot $file.path))
    if (-not @($checkedTargets | Where-Object { $absolute.StartsWith($_ + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase) }).Count) { throw 'Inventory file outside targets' }
    if ((Test-Path -LiteralPath $absolute) -and (Get-FileHash -LiteralPath $absolute -Algorithm SHA256).Hash -ne $file.sha256) { throw 'An obsolete file changed after backup' }
}
# All absolute targets and all backup-bound bytes have now been checked.
foreach ($absolute in $checkedTargets) {
    if (Test-Path -LiteralPath $absolute) {
        Remove-Item -LiteralPath $absolute -Recurse -Force -ErrorAction SilentlyContinue
    }
}
$remaining = @()
foreach ($absolute in $checkedTargets) {
    if (Test-Path -LiteralPath $absolute) {
        $remaining += @(Get-ChildItem -LiteralPath $absolute -Recurse -File -Force | ForEach-Object {
            [ordered]@{ path = [IO.Path]::GetRelativePath($repoRoot, $_.FullName).Replace('\', '/'); bytes = $_.Length }
        })
    }
}
$remainingBytes = [long]0
foreach ($entry in $remaining) { $remainingBytes += [long]$entry.bytes }
$receipt = [ordered]@{
    date = '2026-10-26'
    status = $(if ($remaining.Count) { 'completed_except_in_use_or_unavailable_files' } else { 'completed' })
    attempted_directories = @($plan.targets.directory)
    deleted_directories = @($plan.targets | Where-Object { -not (Test-Path -LiteralPath (Join-Path $repoRoot $_.directory)) } | ForEach-Object { $_.directory })
    remaining_files = $remaining
    removed_files = [long](($plan.targets | Measure-Object -Property files -Sum).Sum) - [long]$remaining.Count
    removed_bytes = [long](($plan.targets | Measure-Object -Property bytes -Sum).Sum) - $remainingBytes
    archived_non_QA_files = $plan.archived_file_count
    discarded_regenerable_QA_files = $plan.discardable_QA_files
    history_archive = $plan.archive
    history_archive_sha256 = $plan.archive_sha256
    submission_package_sha256 = $package.sha256
    v3_and_v4_preserved = $true
    scientific_data_results_models_configs_or_Git_history_deleted = $false
}
$receipt | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $auditDirectory 'cleanup_receipt.json') -Encoding utf8
$receipt | ConvertTo-Json -Depth 6
