[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidateNotNullOrEmpty()]
    [string]$Message,

    [Parameter(Mandatory = $true)]
    [ValidateNotNullOrEmpty()]
    [string[]]$Path,

    [switch]$Yes
)

$ErrorActionPreference = 'Stop'

function Invoke-Git {
    param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Arguments)

    & git @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed with exit code $LASTEXITCODE."
    }
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $repoRoot
try {
    $insideWorkTree = (& git rev-parse --is-inside-work-tree).Trim()
    if ($LASTEXITCODE -ne 0 -or $insideWorkTree -ne 'true') {
        throw 'This script must be run from a Git working tree.'
    }

    # Never combine this focused publish with files staged before the script ran.
    & git diff --cached --quiet
    if ($LASTEXITCODE -eq 1) {
        throw 'The index already contains staged changes. Review, commit, or unstage them before using this script.'
    }
    if ($LASTEXITCODE -ne 0) {
        throw 'Unable to inspect the Git index.'
    }

    $branch = (& git branch --show-current).Trim()
    if ([string]::IsNullOrWhiteSpace($branch)) {
        throw 'Cannot publish from a detached HEAD. Check out a named branch first.'
    }

    Invoke-Git fetch origin --prune
    $remoteRef = "refs/remotes/origin/$branch"
    & git show-ref --verify --quiet $remoteRef
    if ($LASTEXITCODE -eq 0) {
        $counts = (& git rev-list --left-right --count "HEAD...origin/$branch").Trim().Split([char[]]" `t", [System.StringSplitOptions]::RemoveEmptyEntries)
        if ($counts.Count -ne 2) {
            throw "Could not determine ahead/behind status for origin/$branch."
        }
        if ([int]$counts[1] -gt 0) {
            throw "This branch is $($counts[1]) commit(s) behind origin/$branch. Pull or rebase before publishing."
        }
    }
    elseif ($LASTEXITCODE -ne 1) {
        throw "Could not inspect origin/$branch."
    }

    $relativePaths = foreach ($candidate in $Path) {
        if ([IO.Path]::IsPathRooted($candidate)) {
            throw "Use repository-relative paths only: $candidate"
        }
        $fullPath = Join-Path $repoRoot $candidate
        if (-not (Test-Path -LiteralPath $fullPath)) {
            throw "Path does not exist: $candidate"
        }
        $resolvedPath = (Resolve-Path -LiteralPath $fullPath).Path
        $relative = [IO.Path]::GetRelativePath($repoRoot, $resolvedPath).Replace('\', '/')
        if ($relative -eq '..' -or $relative.StartsWith('../')) {
            throw "Path is outside the repository: $candidate"
        }
        & git check-ignore --quiet -- $relative
        if ($LASTEXITCODE -eq 0) {
            throw "Path is ignored by Git and will not be published: $relative"
        }
        if ($LASTEXITCODE -ne 1) {
            throw "Could not determine whether this path is ignored: $relative"
        }
        $relative
    }

    Invoke-Git add -- @relativePaths
    & git diff --cached --quiet
    if ($LASTEXITCODE -eq 0) {
        throw 'The selected paths produced no staged changes.'
    }
    if ($LASTEXITCODE -ne 1) {
        throw 'Unable to inspect staged changes.'
    }

    Invoke-Git diff --cached --check
    Write-Host "`nStaged files:"
    Invoke-Git diff --cached --name-status
    Write-Host "`nStaged summary:"
    Invoke-Git diff --cached --stat

    if (-not $Yes) {
        $confirmation = Read-Host "Type PUBLISH to commit and push branch '$branch'"
        if ($confirmation -cne 'PUBLISH') {
            throw 'Publish cancelled. The selected files remain staged for your review.'
        }
    }

    Invoke-Git commit -m $Message
    Invoke-Git push
    Invoke-Git fetch origin $branch

    $localSha = (& git rev-parse HEAD).Trim()
    $remoteSha = ((& git ls-remote --heads origin $branch) -split "`t")[0].Trim()
    if ($LASTEXITCODE -ne 0 -or $remoteSha -ne $localSha) {
        throw "Push completed but the remote SHA could not be verified for origin/$branch."
    }
    $finalCounts = (& git rev-list --left-right --count "HEAD...origin/$branch").Trim()
    if ($finalCounts -ne '0 0') {
        throw "Push completed but local and remote are not in parity: $finalCounts"
    }

    Write-Host "`nPublished $branch at $localSha (local/remote parity: 0 0)."
}
finally {
    Pop-Location
}
