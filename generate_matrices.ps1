param(
    [int]$J = 2,
    [int]$L = 6,
    [int]$P = 7,
    [int]$M2 = 5,
    [int]$N2 = 30,
    [switch]$NoBuild,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$BuildDir = Join-Path $ScriptDir "build-win"
$Constructor = Join-Path $BuildDir "ldpc_constructor.exe"
$MatrixDir = Join-Path $ScriptDir "matrices"
$SpcLength = $N2 - $M2
if ($SpcLength -lt 2) {
    $SpcLength = 2
}
$Suffix = "_J${J}_L${L}_P${P}_spc${SpcLength}_GF256"

Set-Location $ScriptDir
New-Item -ItemType Directory -Force -Path $MatrixDir | Out-Null

Write-Host ""
Write-Host "Quantum LDPC matrix generation"
Write-Host "Project: $ScriptDir"
Write-Host "Args: J=$J L=$L P=$P m2=$M2 n2=$N2"
Write-Host "Output suffix: $Suffix"
Write-Host ""

if (-not (Test-Path $Constructor)) {
    if ($NoBuild) {
        throw "Missing executable: $Constructor"
    }

    Write-Host "ldpc_constructor.exe was not found. Building it first..."
    cmake -S $ScriptDir -B $BuildDir
    cmake --build $BuildDir --target ldpc_constructor
}

$Command = @($Constructor, "$J", "$L", "$P", "$M2", "$N2")

if ($DryRun) {
    Write-Host "Dry run only. Command would be:"
    Write-Host ($Command -join " ")
    exit 0
}

Write-Host "Running constructor..."
& $Constructor $J $L $P $M2 $N2

if ($LASTEXITCODE -ne 0) {
    throw "Matrix generation failed with exit code $LASTEXITCODE"
}

$ExpectedFiles = @(
    "HX${Suffix}.txt",
    "HZ${Suffix}.txt",
    "HX_Phys${Suffix}.txt",
    "HZ_Phys${Suffix}.txt",
    "H1${Suffix}.txt",
    "H2${Suffix}.txt",
    "H3${Suffix}.txt",
    "G2${Suffix}.txt"
)

Write-Host ""
Write-Host "Generated files:"
foreach ($File in $ExpectedFiles) {
    $Path = Join-Path $MatrixDir $File
    if (Test-Path $Path) {
        Write-Host "  OK  matrices\$File"
    } else {
        Write-Host "  MISS matrices\$File"
    }
}

Write-Host ""
Write-Host "Done. Use this suffix for decoding: $Suffix"
