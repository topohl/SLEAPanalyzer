[CmdletBinding()]
param(
    [string]$ReleaseName = "$(Get-Date -Format 'yyyy-MM-dd')_exp9_all-batches_canonical",
    [string]$BundleRoot = 'S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Behavior\SLEAPanalyzer_v2',
    [string]$AssayOutputRoot = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris\sleap_output_all',
    [string]$StageRoot = 'C:\Users\topohl\Documents\exp9_sleap_staging',
    [string]$HistoricalProjectRoot = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris',
    [string]$EpmCalibrationInput = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris\sleap_input\EPM',
    [string]$NorCalibrationInput = 'C:\Users\topohl\iCloudDrive\Dokumente\Analysis\Behavior\correlate_sleap_boris\sleap_input\NOR',
    [string]$BorisNorRoot = 'S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Raw Data\Behavior\B1\NOR\BORIS',
    [string]$Rscript = 'C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe'
)

$ErrorActionPreference = 'Stop'

if ($ReleaseName -notmatch '^[0-9]{4}-[0-9]{2}-[0-9]{2}_[A-Za-z0-9._-]+$') {
    throw "ReleaseName must start YYYY-MM-DD_ and contain only letters, numbers, dot, underscore or hyphen."
}

$repoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$releaseRoot = Join-Path (Join-Path $BundleRoot 'releases') $ReleaseName

foreach ($required in @(
    $Rscript,
    $AssayOutputRoot,
    $StageRoot,
    $HistoricalProjectRoot,
    $EpmCalibrationInput,
    $NorCalibrationInput,
    $BorisNorRoot
)) {
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
    (Join-Path $releaseRoot 'source_data\validation'),
    (Join-Path $releaseRoot 'source_data\validation\boris_nor'),
    (Join-Path $releaseRoot 'provenance'),
    (Join-Path $releaseRoot 'provenance\configs'),
    (Join-Path $releaseRoot 'provenance\logs'),
    (Join-Path $releaseRoot 'provenance\scripts'),
    (Join-Path $releaseRoot 'provenance\scripts\analysis_core')
)
foreach ($dir in $dirs) {
    New-Item -ItemType Directory -Path $dir -Force | Out-Null
}

$sourceCopies = [System.Collections.Generic.List[object]]::new()
$externalInputs = [System.Collections.Generic.List[object]]::new()
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

function Add-ExternalInput {
    param(
        [Parameter(Mandatory = $true)][string]$Role,
        [Parameter(Mandatory = $true)][string]$Path
    )
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Required external input is missing: $Path"
    }
    $externalInputs.Add([pscustomobject]@{
        role = $Role
        source = $Path
        bytes = (Get-Item -LiteralPath $Path).Length
        sha256 = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
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
    EXP9_SLEAP_SOURCE_DATA_DIR = (Join-Path $releaseRoot 'source_data\validation')
    EXP9_SLEAP_ASSAY_OUTPUT_ROOT = $AssayOutputRoot
    EXP9_SLEAP_STAGE_ROOT = $StageRoot
    EXP9_BORIS_NOR_DIR = $BorisNorRoot
    EXP9_EPM_CALIBRATION_INPUT = $EpmCalibrationInput
    EXP9_NOR_CALIBRATION_INPUT = $NorCalibrationInput
    SLEAPANALYZER_ROOT = $repoRoot
}

try {
    foreach ($name in $environmentUpdates.Keys) {
        $oldEnvironment[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
        [Environment]::SetEnvironmentVariable($name, $environmentUpdates[$name], 'Process')
    }

    $scriptRoot = Join-Path $repoRoot 'validation\exp9_boris\scripts'
    $correlate = Join-Path $scriptRoot '05_correlate.R'
    $noseDipCalibration = Join-Path $scriptRoot '11_recalibrate_nosedips.R'
    $norContactCalibration = Join-Path $scriptRoot '12_calibrate_nor_contact.R'
    $norDetectorComparison = Join-Path $scriptRoot '13_compare_nor_detectors.R'
    $assemble = Join-Path $scriptRoot '09_assemble_all_batches.R'
    $analyse = Join-Path $scriptRoot '10_analyse_all_batches.R'
    $finalise = Join-Path $scriptRoot '11_finalize_output_bundle.R'

    Invoke-RScriptLogged -ScriptPath $correlate -LogName '05_correlate.log'
    Invoke-RScriptLogged -ScriptPath $noseDipCalibration -LogName '11_recalibrate_nosedips.log'
    Invoke-RScriptLogged -ScriptPath $norContactCalibration -LogName '12_calibrate_nor_contact.log'
    Invoke-RScriptLogged -ScriptPath $norDetectorComparison -LogName '13_compare_nor_detectors.log'
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
    foreach ($name in 'epm_all.yaml','nor_all.yaml','socp_all.yaml','oft_all.yaml','epm_b1.yaml','nor_b1.yaml') {
        Copy-RecordedFile `
            -Source (Join-Path $configDir $name) `
            -Destination (Join-Path $releaseRoot "provenance\configs\$name")
    }
    foreach ($name in '00_theme.R','05_correlate.R','08_stage_all_batches.R','09_assemble_all_batches.R','10_analyse_all_batches.R','11_finalize_output_bundle.R','11_recalibrate_nosedips.R','12_calibrate_nor_contact.R','13_compare_nor_detectors.R') {
        Copy-RecordedFile `
            -Source (Join-Path $repoRoot "validation\exp9_boris\scripts\$name") `
            -Destination (Join-Path $releaseRoot "provenance\scripts\$name")
    }
    Copy-RecordedFile `
        -Source $PSCommandPath `
        -Destination (Join-Path $releaseRoot 'provenance\scripts\build_exp9_behavior_bundle.ps1')

    $validationRoot = Join-Path $repoRoot 'validation\exp9_boris'
    Copy-RecordedFile `
        -Source (Join-Path $validationRoot 'metadata\novelLoc.txt') `
        -Destination (Join-Path $releaseRoot 'source_data\validation\novelLoc.txt')
    Get-ChildItem -LiteralPath $BorisNorRoot -File -Filter '*_nov.tsv' | Sort-Object Name | ForEach-Object {
        Copy-RecordedFile `
            -Source $_.FullName `
            -Destination (Join-Path $releaseRoot "source_data\validation\boris_nor\$($_.Name)")
    }

    $analysisCode = @(
        '02_SLEAPanalzyer\DLCAnalyzer_Functions_final.R',
        '02_SLEAPanalzyer\Behavioral_Metrics_Phase1.R',
        '02_SLEAPanalzyer\core\events.R',
        '02_SLEAPanalzyer\core\geometry.R',
        '02_SLEAPanalzyer\core\interpolation.R',
        '02_SLEAPanalzyer\core\io.R',
        '02_SLEAPanalzyer\core\validation.R'
    )
    foreach ($relative in $analysisCode) {
        Copy-RecordedFile `
            -Source (Join-Path $repoRoot $relative) `
            -Destination (Join-Path $releaseRoot "provenance\scripts\analysis_core\$(Split-Path -Leaf $relative)")
    }
    Copy-RecordedFile `
        -Source (Join-Path $repoRoot '02_SLEAPanalzyer\EPM_zoneinfo.csv') `
        -Destination (Join-Path $releaseRoot 'source_data\validation\EPM_zoneinfo.csv')

    Add-ExternalInput `
        -Role 'Method validation manual source table' `
        -Path (Join-Path $validationRoot 'enriched\analysis_ready_wide.tsv')
    Add-ExternalInput `
        -Role 'Method validation SLEAP source table' `
        -Path (Join-Path $validationRoot 'enriched\sleap_wide.tsv')
    Get-ChildItem -LiteralPath $EpmCalibrationInput -File -Filter '*.csv' | Sort-Object Name | ForEach-Object {
        Add-ExternalInput -Role 'EPM calibration coordinates' -Path $_.FullName
    }
    Get-ChildItem -LiteralPath $NorCalibrationInput -File -Filter '*.csv' | Sort-Object Name | ForEach-Object {
        Add-ExternalInput -Role 'NOR calibration coordinates' -Path $_.FullName
    }
    $externalInputs | Export-Csv `
        -LiteralPath (Join-Path $releaseRoot 'provenance\external_input_manifest.tsv') `
        -Delimiter "`t" -NoTypeInformation -Encoding utf8

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
        "epm_calibration_input=$EpmCalibrationInput",
        "nor_calibration_input=$NorCalibrationInput",
        "boris_nor_root=$BorisNorRoot",
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
