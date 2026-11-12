"""`launch_detached` の中継の代わりに眠るだけのモジュール（試験専用）。

`tests/test_process_control.py` が `process_control._RELAY_MODULE` をこれに差し替え、
中継が上限までに終わらないときの経路（`DETACH_RELAY_TIMEOUT`）を見る。
後片付けに失敗しても 60 秒で自分から終わる。
"""
import time

time.sleep(60)
