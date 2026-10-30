# 再解析 v8: Console の試し（Task 10）の記録

2026-10-30〜31。MTBLS6965 の R17_1-3・P17_1-3（陽陰 6 ファイルずつ）を `library_mode=msp_then_lbm` で回した。
MS-DIAL Console は MsdialCoreTestApp（Debug net48、45a531c4d）。作業場所は `datasets/reanalysis/_v8_trial/`。

## 1 回目（修正前、31c4b6f）

- 未署名の `MsdialGcMsApi.dll` が Smart App Control に拒否され、2 本とも failed になった（CodeIntegrity 3033/3077）。
  PC の再起動の後、neg は通った。pos は再起動の直後にもう一度拒否され、次の再開で通った。
- 2 本とも completed になった。2 本を同時に回し、約 33 分かかった。
  MSDIALCUI の WorkingSet は、1 本の最大が 25.8 GB、2 本の合計の最大が 27.1 GB だった。
- 結果は期待と違った。研究室 MSP の候補が足切りを越えるのに LBM が代表になる特徴が、26 件あった（pos 12、neg 14）。
  例は TG 54:3 と Glyceryl trioleate。
- 原因: Console は `Weighted/Simple/Reverse dot product cutoff for MSP-based annotation` の値を**二乗値**として読む。
  そのため、メソッドの 0.6/0.6/0.8 は実効で 0.775/0.775/0.894 になる。
  LBM には足切りのキーが無いので既定の 0.6/0.6/0.8 が使われ、MSP だけが厳しくなっていた。
  MSP だけの 1 本（R17_3 pos）でこの規則を確かめた。.mdpeak の `MS/MS matched` は、13 ファイルで食い違いが 0 だった。
- 修正 de70012: `msp_then_lbm` は、MSP と LBM の dot product の足切りを `Square root of ... cutoff for {MSP,LBM}-based annotation` のキーで同じ値に書き、二乗のキーは空にする。

## 回し直し（de70012）

| 確かめること | pos | neg |
|---|---|---|
| 完走 | completed | completed |
| アライメントの基準（`alignment_reference`） | file_id 2 = R17_3（P17-R-33） | file_id 2 = R17_3（P17-R-33） |
| effective-method | 5 キーと平方根の足切り 6 キー（0.6/0.6/0.8）。二乗のキー 3 つは空 | 同じ |
| mzTab の SME（MSP / LBM） | 835 / 528（1 回目 804 / 559） | 602 / 571（1 回目 563 / 610） |
| 代表が MSP で matched | 553（1 回目 448） | 356（1 回目 251） |
| 代表が LBM で matched | 261（1 回目 292） | 223（1 回目 262） |
| LBM が代表なのに、MSP の候補が足切りを越える特徴 | 0（1 回目 12） | 1（1 回目 14） |

effective-method の 5 キーは、`Msp file path` が空、`MSP annotator settings file path`、`Lbm file path`、LBM の RT のキー 2 つが False。

- **両方の候補があるとき、rank 1 が MSP になる例**: pos の Alignment ID 1284（m/z 902.822、RT 11.99、[M+NH4]+）。
  1 回目の代表は LBM の `TG 54:3|TG 18:0_18:1_18:2` だった。
  回し直しでは研究室 MSP の `Glyceryl trioleate`（matched、wdp 0.666・sdp 0.712・rdp 0.863）になった。
- 残る 1 件は neg の 1108（`PC 35:1|PC 17:0_18:1`）。こちらで採点し直すと、MSP の PC 35:1 は wdp 0.748 で足切りを越える。
  ただし、この採点は .mdalign の代表スペクトルに対するもので、MS-DIAL の判定（代表ピークのスペクトル）とは少しずれ得る。
  1 回目に採点が一致したのは 593 件中 592 件だった。この 1 件はそのずれの範囲とみなし、追っていない。
- メモリ: 回し直しの最大は記録できていない（見張りのシェルがメモリ不足で止められた）。設定はキーの値しか変わらないので、1 回目の値を使う。

## 残り

- 過去の accession（MTBLS3233、MTBLS6039、MTBLS6965、ST002901、ST003841）の effective-method にも、同じ二乗のキー（0.6/0.6/0.8）がある。
  そのため、研究室 MSP の同定は意図より厳しかった可能性がある。どう扱うかは別に決める。
