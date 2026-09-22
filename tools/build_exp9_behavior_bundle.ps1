[CmdletBinding()]
param(
    [string]$ReleaseName = "$(Get-Date -Format 'yyyy-MM-dd')_exp9_all-batches_canonical",
    [string]$BundleRoot = 'S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Behavior\SLEAPanalyzer_v2',
    [string]$AssayOutputRoot = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris\sleap_output_all',
    [string]$StageRoot = 'C:\Users\topohl\Documents\exp9_sleap_staging',
    [string]$HistoricalProjectRoot = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris',
    [string]$Rscript = 'C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe'
)

$ErrorActionPreference = 'Stop'

if ($ReleaseName -notmatch '^[0-9]{4}-[0-9]{2}-[0-9]{2}_[A-Za-z0-9._-]+$') {
    throw "ReleaseName must start YYYY-MM-DD_ and contain only letters, numbers, dot, underscore or hyphen."
}

$repoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$releaseRoot = Join-Path (Join-Path $BundleRoot 'releases') $ReleaseName

foreach ($required in @($Rscript, $AssayOutputRoot, $StageRoot, $HistoricalProjectRoot)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Required path does not exist: $required"
    }
}
if (Test-Path -LiteralPath $releaseRoot) {
    throw "Release already exists and will not be overwritten: $releaseRoot"
}

$dirs = @(
    $releaseRoot,
    (Join-Path $releaseRoot 'data'),
    (Join-Path $releaseRoot 'metadata'),
    (Join-Path $releaseRoot 'assay_summaries'),
    (Join-Path $releaseRoot 'statistics'),
    (Join-Path $releaseRoot 'figures'),
    (Join-Path $releaseRoot 'qc'),
    (Join-Path $releaseRoot 'provenance'),
    (Join-Path $releaseRoot 'provenance\configs'),
    (Join-Path $releaseRoot 'provenance\logs'),
    (Join-Path $releaseRoot 'provenance\scripts')
)
foreach ($dir in $dirs) {
    New-Item -ItemType Directory -Path $dir -Force | Out-Null
}

$sourceCopies = [System.Collections.Generic.List[object]]::new()
function Copy-RecordedFile {
    param(
        [Parameter(Mandatory = $true)][string]$Source,
        [Parameter(Mandatory = $true)][string]$Destination
    )
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        throw "Required bundle source is missing: $Source"
    }
    $parent = Split-Path -Parent $Destination
    New-Item -ItemType Directory -Path $parent -Force | Out-Null
    Copy-Item -LiteralPath $Source -Destination $Destination
    $sourceHash = (Get-FileHash -LiteralPath $Source -Algorithm SHA256).Hash.ToLowerInvariant()
    $destinationHash = (Get-FileHash -LiteralPath $Destination -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($sourceHash -ne $destinationHash) {
        throw "Hash mismatch after copying $Source"
    }
    $sourceCopies.Add([pscustomobject]@{
        source = $Source
        bundled_path = $Destination.Substring($releaseRoot.Length + 1)
        bytes = (Get-Item -LiteralPath $Destination).Length
        sha256 = $destinationHash
    })
}

function Invoke-RScriptLogged {
    param(
        [Parameter(Mandatory = $true)][string]$ScriptPath,
        [Parameter(Mandatory = $true)][string]$LogName
    )
    $logPath = Join-Path $releaseRoot "provenance\logs\$LogName"
    & $Rscript $ScriptPath *>&1 | Tee-Object -FilePath $logPath
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) {
        throw "$(Split-Path -Leaf $ScriptPath) failed with exit code $exitCode"
    }
}

$oldEnvironment = @{}
$environmentUpdates = @{
    EXP9_SLEAP_RUN_ROOT = $releaseRoot
    EXP9_SLEAP_DATA_DIR = (Join-Path $releaseRoot 'data')
    EXP9_SLEAP_METADATA_DIR = (Join-Path $releaseRoot 'metadata')
    EXP9_SLEAP_RESULTS_DIR = (Join-Path $releaseRoot 'statistics')
    EXP9_SLEAP_FIGURES_DIR = (Join-Path $releaseRoot 'figures')
    EXP9_SLEAP_ASSAY_OUTPUT_ROOT = $AssayOutputRoot
    EXP9_SLEAP_STAGE_ROOT = $StageRoot
}

try {
    foreach ($name in $environmentUpdates.Keys) {
        $oldEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
        [Environment]::SetEnvironmentVariable($name, $environmentUpdates[$name], 'Process')
    }

    $assemble = Join-Path $repoRoot 'validation\exp9_boris\scripts\09_assemble_all_batches.R'
    $analyse = Join-Path $repoRoot 'validation\exp9_boris\scripts\10_analyse_all_batches.R'
    $finalise = Join-Path $repoRoot 'validation\exp9_boris\scripts\11_finalize_output_bundle.R'

    Invoke-RScriptLogged -ScriptPath $assemble -LogName '09_assemble_all_batches.log'
    Invoke-RScriptLogged -ScriptPath $analyse -LogName '10_analyse_all_batches.log'
    Invoke-RScriptLogged -ScriptPath $finalise -LogName '11_finalize_output_bundle.log'

    $batches = 'B1','B2','B3','B4','B5','B6'
    foreach ($batch in $batches) {
        foreach ($name in 'Report.csv','tracking_qc.csv','run_manifest.yaml') {
            $source = Join-Path $AssayOutputRoot "EPM\$batch\$name"
            $destination = Join-Path $releaseRoot "assay_summaries\EPM\$batch\$name"
            Copy-RecordedFile -Source $source -Destination $destination
        }
        foreach ($name in 'combined_output.csv','tracking_qc.csv','run_manifest.yaml') {
            $source = Join-Path $AssayOutputRoot "NOR\$batch\$name"
            $destination = Join-Path $releaseRoot "assay_summaries\NOR\$batch\$name"
            Copy-RecordedFile -Source $source -Destination $destination
        }
        foreach ($phase in 'S1','S2') {
            foreach ($name in 'combined_output.csv','tracking_qc.csv','run_manifest.yaml') {
                $source = Join-Path $AssayOutputRoot "SocP\$batch\$phase\$name"
                $destination = Join-Path $releaseRoot "assay_summaries\SocP\$batch\$phase\$name"
                Copy-RecordedFile -Source $source -Destination $destination
            }
        }
        foreach ($name in 'combined_summary.csv','combined_qc.csv','combined_bins.csv','combined_enhanced_metrics.csv','combined_first_last_epoch_metrics.csv') {
            $source = Join-Path $StageRoot "OFT\$batch\output_v1.2.0\$name"
            $destination = Join-Path $releaseRoot "assay_summaries\OFT\$batch\$name"
            Copy-RecordedFile -Source $source -Destination $destination
        }
    }

    $qcSources = 'staging_report.csv','batch_audit_files.csv','batch_audit_geometry_epm.csv','batch_audit_geometry_rect.csv'
    foreach ($name in $qcSources) {
        Copy-RecordedFile `
            -Source (Join-Path $HistoricalProjectRoot "results\$name") `
            -Destination (Join-Path $releaseRoot "qc\source_$name")
    }

    $configDir = Join-Path $repoRoot 'validation\exp9_boris\config'
    foreach ($name in 'epm_all.yaml','nor_all.yaml','socp_all.yaml','oft_all.yaml') {
        Copy-RecordedFile `
            -Source (Join-Path $configDir $name) `
            -Destination (Join-Path $releaseRoot "provenance\configs\$name")
    }
    foreach ($name in '00_theme.R','08_stage_all_batches.R','09_assemble_all_batches.R','10_analyse_all_batches.R','11_finalize_output_bundle.R') {
        Copy-RecordedFile `
            -Source (Join-Path $repoRoot "validation\exp9_boris\scripts\$name") `
            -Destination (Join-Path $releaseRoot "provenance\scripts\$name")
    }
    Copy-RecordedFile `
        -Source $PSCommandPath `
        -Destination (Join-Path $releaseRoot 'provenance\scripts\build_exp9_behavior_bundle.ps1')

    $sourceCopies | Export-Csv `
        -LiteralPath (Join-Path $releaseRoot 'provenance\source_copy_manifest.tsv') `
        -Delimiter "`t" -NoTypeInformation -Encoding utf8

    $gitLines = @(
        "repository=$repoRoot",
        "commit=$(& git -C $repoRoot rev-parse HEAD)",
        "branch=$(& git -C $repoRoot branch --show-current)",
        'status:',
        (& git -C $repoRoot status --short)
    )
    [System.IO.File]::WriteAllLines(
        (Join-Path $releaseRoot 'provenance\git_state.txt'),
        [string[]]$gitLines
    )
    $parameterLines = @(
        "release_name=$ReleaseName",
        "release_root=$releaseRoot",
        "assay_output_root=$AssayOutputRoot",
        "stage_root=$StageRoot",
        "historical_project_root=$HistoricalProjectRoot",
        "rscript=$Rscript"
    )
    [System.IO.File]::WriteAllLines(
        (Join-Path $releaseRoot 'provenance\build_parameters.txt'),
        [string[]]$parameterLines
    )
    & $Rscript -e "writeLines(capture.output(sessionInfo()), '$($releaseRoot.Replace('\','/'))/provenance/r_session_info.txt')"
    if ($LASTEXITCODE -ne 0) { throw "Could not capture R session information" }

    $manifestRows = Get-ChildItem -LiteralPath $releaseRoot -Recurse -File |
        Where-Object { $_.Name -ne 'bundle_manifest.tsv' } |
        Sort-Object FullName |
        ForEach-Object {
            [pscustomobject]@{
                path = $_.FullName.Substring($releaseRoot.Length + 1)
                bytes = $_.Length
                sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
            }
        }
    $manifestRows | Export-Csv `
        -LiteralPath (Join-Path $releaseRoot 'provenance\bundle_manifest.tsv') `
        -Delimiter "`t" -NoTypeInformation -Encoding utf8

    Write-Output "Bundle created: $releaseRoot"
    Write-Output "Files in manifest: $($manifestRows.Count)"
}
catch {
    [System.IO.File]::WriteAllText(
        (Join-Path $releaseRoot 'INCOMPLETE.txt'),
        "Bundle construction failed: $($_.Exception.Message)"
    )
    throw
}
finally {
    foreach ($name in $environmentUpdates.Keys) {
        [Environment]::SetEnvironmentVariable($name, $oldEnvironment[$name], 'Process')
    }
}
