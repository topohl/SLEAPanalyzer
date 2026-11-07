# Exp9 quality checks in one go: the video QC manifest (32), the video checks (tools/video_qc), the plausibility screen
# (23) when it is older than the staged tracking, the light checks (31), then a review copy on the lab share. Report only:
# nothing here changes a score.
#
#   powershell -File tools\run_exp9_qc.ps1                    # frames and objects, about 30 min in all
#   powershell -File tools\run_exp9_qc.ps1 -Checks sync       # add a slower check; the others' results are kept
#   powershell -File tools\run_exp9_qc.ps1 -SkipVideo         # light checks only (about 10-20 min)
#   powershell -File tools\run_exp9_qc.ps1 -WhatIf            # name the steps it would run, and why, without running them
#
# build_exp9_behavior_bundle.ps1 runs it at the end of every release build (-Release, -ReleaseQcDir).
#
# 31 joins results\plausibility_summary.csv of 23_plausibility_screen.R. When that file is missing or older than a staged
# tracking file, 23 runs first (about an hour more); -SkipPlausibility leaves it, and the warnings at the end say so.
#
# At the end it lists the flags still to review: those without a decision in
# validation\exp9_boris\config\light_check_decisions.csv, or on a tracking file changed since. Each such recording gets
# a review sheet (tools\video_qc\review.py) in video_qc\stills\review. After looking at one, add its row to the decisions
# file (outcome keep, correct, code_fixed or open; finding; evidence) and run again. The video checks are incremental:
# a run checks only recordings whose tracking, video or settings changed, and continues a run that stopped.
[CmdletBinding(SupportsShouldProcess)]
param(
    [string]$Checks = 'frames,objects',
    [switch]$SkipVideo,
    [switch]$SkipPlausibility,      # do not rerun 23 even when its summary is older than the tracking
    [switch]$NoCopy,
    [string]$Release = '',          # the release 31 reads (default: the one named in 31)
    [string]$StageRoot = '',        # the staged tracking (default: EXP9_STAGING, else the one named in 23, 31 and 32)
    [string]$ReleaseQcDir = '',     # also copy the results here (a release's qc folder)
    [string]$ShareDir = 'S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Behavior\SLEAP\qc\light_checks',
    [string]$Rscript = 'C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe',
    [string]$Python = 'C:\Users\topohl\AppData\Local\miniconda3\python.exe'
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$study = Join-Path $repo 'validation\exp9_boris'
$results = Join-Path $study 'results'
$scripts = Join-Path $study 'scripts'
$warnings = [System.Collections.Generic.List[string]]::new()   # repeated at the end

function Invoke-Step([string]$What, [scriptblock]$Run, [string]$Why = '') {
    Write-Host "`n== $What" -ForegroundColor Cyan
    if ($Why) { Write-Host "   $Why" }
    if (-not $PSCmdlet.ShouldProcess($What, 'run')) { return }
    $t = Get-Date
    & $Run
    if ($LASTEXITCODE -ne 0) { throw "$What failed (exit $LASTEXITCODE)" }
    Write-Host ("   {0:n0} min" -f ((Get-Date) - $t).TotalMinutes)
}

function Copy-Results([string]$To) {
    Write-Host "`n== copy to $To" -ForegroundColor Cyan
    if (-not $PSCmdlet.ShouldProcess($To, 'copy the results')) { return }
    $vq = Join-Path $To 'video_qc'
    New-Item -ItemType Directory -Force -Path $vq | Out-Null
    Copy-Item (Join-Path $results 'light_checks_*.csv') $To -Force
    Copy-Item (Join-Path $study 'config\light_check_decisions.csv') $To -Force
    Copy-Item (Join-Path $results 'video_qc\*.csv') $vq -Force
    Copy-Item (Join-Path $results 'video_qc_manifest.csv') $vq -Force
    $stills = Join-Path $results 'video_qc\stills'
    if (Test-Path -LiteralPath $stills) {
        robocopy $stills (Join-Path $vq 'stills') /MIR /NJH /NJS /NFL /NDL /NP | Out-Null   # stills of this run only
        if ($LASTEXITCODE -ge 8) { throw "copying the stills failed (robocopy $LASTEXITCODE)" }
        $global:LASTEXITCODE = 0
    }
}

function Get-StagedTracking([string]$Root) {
    # the tracking files 23, 31 and 32 list: csv files directly in these folders, not in the output folders below them
    $dirs = foreach ($b in 1..6) {
        foreach ($p in 'HAB', 'S1', 'S2') { "SocP\B$b\$p" }
        "NOR\B$b", "OFT\B$b", "EPM\B$b"
        if ($b -gt 1) { "NOR_HAB\B$b" }
    }
    foreach ($d in $dirs) {
        $p = Join-Path $Root $d
        if (Test-Path -LiteralPath $p -PathType Container) {
            Get-ChildItem -LiteralPath $p -File | Where-Object { $_.Name -cmatch '[.]csv$' }
        }
    }
}

function Get-PlausibilityStale([string]$Summary, [object[]]$Tracking) {
    # why 23 should run before 31 (its summary missing, or older than a staged tracking file); nothing when up to date
    if (-not (Test-Path -LiteralPath $Summary)) { return "no $Summary" }
    $written = (Get-Item -LiteralPath $Summary).LastWriteTimeUtc
    $newer = @($Tracking | Where-Object { $_.LastWriteTimeUtc -gt $written } | Sort-Object LastWriteTimeUtc -Descending)
    if ($newer.Count) {
        '{0} of {1} staged tracking files are newer than plausibility_summary.csv ({2:yyyy-MM-dd HH:mm}), the newest {3} ({4:yyyy-MM-dd HH:mm})' -f `
            $newer.Count, $Tracking.Count, $written.ToLocalTime(), $newer[0].FullName, $newer[0].LastWriteTime
    }
}

# the staging root the R scripts read: -StageRoot, else EXP9_STAGING, else their default (as written in 23)
$stage = if ($StageRoot) { $StageRoot } elseif ($env:EXP9_STAGING) { $env:EXP9_STAGING } else {
    $m = Select-String -LiteralPath (Join-Path $scripts '23_plausibility_screen.R') -Pattern 'Sys\.getenv\("EXP9_STAGING", "([^"]+)"\)' |
        Select-Object -First 1
    if (-not $m) { throw "no EXP9_STAGING default found in 23_plausibility_screen.R: give -StageRoot" }
    $m.Matches[0].Groups[1].Value
}
$tracking = @(Get-StagedTracking $stage)
if (-not $tracking.Count) { throw "no staged tracking files under $stage" }
Write-Host "staged tracking: $($tracking.Count) files under $stage"

$oldEnv = @{ EXP9_RELEASE = $env:EXP9_RELEASE; EXP9_STAGING = $env:EXP9_STAGING }
try {
    if ($Release) { $env:EXP9_RELEASE = $Release }
    if ($StageRoot) { $env:EXP9_STAGING = $StageRoot }

    Invoke-Step 'video QC manifest (32)' { & $Rscript (Join-Path $scripts '32_video_qc_manifest.R') }
    if (-not $SkipVideo) {
        Invoke-Step "video checks: $Checks" {
            & $Python -I (Join-Path $repo 'tools\video_qc') (Join-Path $results 'video_qc_manifest.csv') `
                (Join-Path $results 'video_qc') --checks $Checks --stills
        }
    }

    $stale = Get-PlausibilityStale (Join-Path $results 'plausibility_summary.csv') $tracking
    if ($stale -and $SkipPlausibility) {
        Write-Host "`n== plausibility screen (23) not rerun (-SkipPlausibility)" -ForegroundColor Cyan
        Write-Host "   $stale"
        $warnings.Add("plausibility screen (23) not rerun (-SkipPlausibility): $stale; the jump, skeleton, outside_arena and raised values of 31 may be stale")
    }
    elseif ($stale) {
        # a failure leaves the summary as it was: 31 still runs and warns that those values may be stale
        try {
            Invoke-Step 'plausibility screen (23)' { & $Rscript (Join-Path $scripts '23_plausibility_screen.R') } `
                -Why "$stale; rerunning it (about an hour)"
        }
        catch {
            Write-Warning $_.Exception.Message
            $warnings.Add("$($_.Exception.Message): 31 used the plausibility summary from before ($stale)")
        }
    }
    Invoke-Step 'light checks (31)' { & $Rscript (Join-Path $scripts '31_light_checks.R') }

    # a review sheet for every recording with flags to review (video_qc\stills\review, copied with the other stills)
    try {
        Invoke-Step 'review sheets' {
            & $Python -I (Join-Path $repo 'tools\video_qc\review.py') (Join-Path $results 'video_qc_manifest.csv') `
                (Join-Path $results 'light_checks_recordings.csv') (Join-Path $results 'video_qc\stills\review') `
                --id-columns 'assay phase batch Code' --draw to_review --pairs (Join-Path $results 'light_checks_duplicates.csv')
        }
    }
    catch {
        Write-Warning $_.Exception.Message
        $warnings.Add("$($_.Exception.Message): no review sheets this time")
    }

    if ($ReleaseQcDir) { Copy-Results $ReleaseQcDir }
    if (-not $NoCopy) { Copy-Results $ShareDir }

    if (-not $WhatIfPreference) {
        $todo = Import-Csv (Join-Path $results 'light_checks_recordings.csv') | Where-Object { $_.to_review }
        Write-Host "`n== flags to review: $(@($todo).Count) recordings" -ForegroundColor Cyan
        if ($todo) { $todo | Format-Table assay, phase, batch, Code, to_review -AutoSize -Wrap | Out-String | Write-Host }
        else { Write-Host 'none: every flag has a decision' }
    }
    if ($warnings.Count) {
        Write-Host "`n== warnings" -ForegroundColor Yellow
        foreach ($w in $warnings) { Write-Host "   $w" }
    }
}
finally {
    # [NullString]::Value removes a variable that was not set; a plain $null would leave it set to "" (PowerShell 7)
    foreach ($k in $oldEnv.Keys) {
        [Environment]::SetEnvironmentVariable($k, $(if ($null -eq $oldEnv[$k]) { [NullString]::Value } else { $oldEnv[$k] }), 'Process')
    }
}
