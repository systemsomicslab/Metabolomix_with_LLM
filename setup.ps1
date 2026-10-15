<#
.SYNOPSIS
  ms-data-parser の Python 環境（.venv）と MCP 設定を 1 回で揃える。

.DESCRIPTION
  リポジトリ直下で実行する:  ./setup.ps1        （利用者）
                             ./setup.ps1 -Dev   （開発者。pytest も入る）
  uv があれば uv.lock どおりに入れる。無ければ Python 3.13 以上で標準の venv を作り、
  requirements.lock.txt（uv.lock の書き出し）を pip で入れる。最後に
  scripts/setup_finish.py が import の自己確認と .mcp.json / lipidmix.local.toml の生成を行う。
  R（パスウェイ投影）は含まない: 使うときだけ Rscript r/setup.R。

.PARAMETER Dev
  開発用の依存（pytest）も入れる。

.PARAMETER NoUv
  uv があっても使わず、標準の venv + pip で作る（検証用）。
#>
param([switch]$Dev, [switch]$NoUv)
# 'Stop' にしない: Windows PowerShell 5.1 は native コマンドの stderr 出力を終了エラーにして
# しまい、$LASTEXITCODE の確認や Fail の 1 行に届かないため。失敗は終了コードで判定する。
$ErrorActionPreference = 'Continue'
Set-Location -LiteralPath $PSScriptRoot -ErrorAction Stop

function Fail([string]$message) {
    Write-Host "[setup] 失敗: $message" -ForegroundColor Red
    exit 1
}

function Find-Python {
    # 3.13 以上の Python を探す。見つかれば起動に使う引数列を返す。
    $candidates = @(@('py', '-3.14'), @('py', '-3.13'), @('python'), @('python3'))
    foreach ($candidate in $candidates) {
        if (-not (Get-Command $candidate[0] -ErrorAction SilentlyContinue)) { continue }
        $exe = $candidate[0]
        $rest = @($candidate | Select-Object -Skip 1)
        & $exe @rest -c "import sys; sys.exit(0 if sys.version_info >= (3, 13) else 1)" 2>$null
        if ($LASTEXITCODE -eq 0) { return ,$candidate }
    }
    return $null
}

$uv = if ($NoUv) { $null } else { Get-Command uv -ErrorAction SilentlyContinue }
if ($uv) {
    Write-Host '[setup] uv で .venv を作ります（uv.lock どおり）。'
    $group = if ($Dev) { @('--group', 'dev') } else { @('--no-dev') }
    & uv sync --locked @group
    if ($LASTEXITCODE -ne 0) { Fail 'uv sync が失敗しました。上の出力を確認してください。' }
} else {
    $python = Find-Python
    if (-not $python) {
        Fail ('Python 3.13 以上も uv も見つかりません。どちらかを入れてから再実行してください。' +
              ' uv: https://docs.astral.sh/uv/getting-started/installation/ ' +
              ' Python: https://www.python.org/downloads/')
    }
    Write-Host "[setup] uv が無いので標準の venv で作ります（$($python -join ' ')）。"
    $exe = $python[0]
    $rest = @($python | Select-Object -Skip 1)
    & $exe @rest -m venv .venv
    if ($LASTEXITCODE -ne 0) { Fail 'python -m venv が失敗しました。' }
    $lock = if ($Dev) { 'requirements-dev.lock.txt' } else { 'requirements.lock.txt' }
    & .\.venv\Scripts\python.exe -m pip install --disable-pip-version-check -r $lock
    if ($LASTEXITCODE -ne 0) { Fail "pip install -r $lock が失敗しました。上の出力を確認してください。" }
}

& .\.venv\Scripts\python.exe scripts\setup_finish.py
if ($LASTEXITCODE -ne 0) { Fail 'セットアップの確認で問題が見つかりました（上の出力）。' }
