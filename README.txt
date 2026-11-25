統合分析セット（架空データ）

1. ZIPを新しいフォルダに展開する。
2. RStudioで02_integrated_analysis.Rを開く。
3. Session > Set Working Directory > To Source File Location を選ぶ。
4. 初回のみ install.packages(c("survival", "ordinal", "nnet")) をConsoleで実行。
5. Sourceを押す。結果はresults_Rに保存される。

既存結果を上書きしないよう、結果フォルダが空でない場合は停止する。
再実行時は結果フォルダを別名で保存するか、新規フォルダを使用する。
01_generate_data.pyは架空データ生成用の原本。既にCSVがあるので通常は実行不要。
元の02〜12の処理順序と推定式を維持し、見出し、パッケージ読み込み、保存処理、安全確認を追加した。
R実行環境がないため動作および元結果との数値一致は未検証。CLMM推論は暫定。
