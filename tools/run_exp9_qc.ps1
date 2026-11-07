# Exp9 quality checks in one go: the video QC manifest (32), the video checks (tools/video_qc), the light checks (31),
# then a review copy on the lab share. Report only: nothing here changes a score.
#
#   powershell -File tools\run_exp9_qc.ps1                    # frames and objects, about 30 min in all
#   powershell -File tools\run_exp9_qc.ps1 -Checks sync       # add a slower check; the others' results are kept
#   powershell -File tools\run_exp9_qc.ps1 -SkipVideo         # light checks only (about 10-20 min)
#
# build_exp9_behavior_bundle.ps1 runs it at the end of every release build (-Release, -ReleaseQcDir).
#
# At the end it lists the flags still to review: those without a decision in
# validation\exp9_boris\config\light_check_decisions.csv, or on a tracking file changed since. After looking at one,
# add its row there (outcome keep, correct, code_fixed or open; finding; evidence) and run again.
[CmdletBinding()]
param(
    [string]$Checks = 'frames,objects',
    [switch]$SkipVideo,
    [switch]$NoCopy,
    [string]$Release = '',          # the release 31 reads (default: the one named in 31)
    [string]$StageRoot = '',        # the staged tracking (default: the one named in 31 and 32)
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

function Invoke-Step([string]$What, [scriptblock]$Run) {
    Write-Host "`n== $What" -ForegroundColor Cyan
    $t = Get-Date
    & $Run
    if ($LASTEXITCODE -ne 0) { throw "$What failed (exit $LASTEXITCODE)" }
    Write-Host ("   {0:n0} min" -f ((Get-Date) - $t).TotalMinutes)
}

function Copy-Results([string]$To) {
    Write-Host "`n== copy to $To" -ForegroundColor Cyan
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
    Invoke-Step 'light checks (31)' { & $Rscript (Join-Path $scripts '31_light_checks.R') }

    if ($ReleaseQcDir) { Copy-Results $ReleaseQcDir }
    if (-not $NoCopy) { Copy-Results $ShareDir }

    $todo = Import-Csv (Join-Path $results 'light_checks_recordings.csv') | Where-Object { $_.to_review }
    Write-Host "`n== flags to review: $(@($todo).Count) recordings" -ForegroundColor Cyan
    if ($todo) { $todo | Format-Table assay, phase, batch, Code, to_review -AutoSize -Wrap | Out-String | Write-Host }
    else { Write-Host 'none: every flag has a decision' }
}
finally {
    foreach ($k in $oldEnv.Keys) { [Environment]::SetEnvironmentVariable($k, $oldEnv[$k], 'Process') }
}
