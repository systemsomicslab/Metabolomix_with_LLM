import os
import io
import msgpack
import lz4.block
from metabolomix.arf.reader import deserialize_lz4_packed_msgpack

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
path = 'data/AlignmentResult_2026_04_07_14_07_25.arf2'
with open(path, 'rb') as f:
    data = f.read()
all_data = deserialize_lz4_packed_msgpack(data)
spot = all_data[0][3][0]
print('spot len', len(spot))
for idx in [4, 10, 41, 42, 53, 56]:
    v = spot[idx]
    print('idx', idx, 'type', type(v).__name__, 'len', len(v) if hasattr(v, '__len__') else None)
    if isinstance(v, list) and len(v) > 0:
        print('  first element type', type(v[0]).__name__, 'len', len(v[0]) if hasattr(v[0], '__len__') else None)
        print('  repr', repr(v[0])[:300])
    print('---')
