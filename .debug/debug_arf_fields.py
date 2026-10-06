import os
import io
import msgpack
import lz4.block
from metabolomix.arf.reader import deserialize_lz4_packed_msgpack

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
path = r'data\AlignmentResult_2026_04_07_14_07_25_PeakProperties.arf'
with open(path, 'rb') as f:
    data = f.read()
all_data = deserialize_lz4_packed_msgpack(data)
print('top count', len(all_data))
for i, d in enumerate(all_data):
    print('--- top', i, 'type', type(d).__name__, 'len', len(d))
    if i == 0:
        print('first 10 types', [type(x).__name__ for x in d[:10]])
        print('d0[0]', d[0])
        print('d0[1]', d[1])
        print('d0[2] type', type(d[2]).__name__, 'len', len(d[2]))
        print('d0[2][0] type', type(d[2][0]).__name__, 'len', len(d[2][0]) if hasattr(d[2][0], '__len__') else None)
        print('d0[2][0] repr', repr(d[2][0])[:300])
        print('d0[3] type', type(d[3]).__name__, 'len', len(d[3]))
        print('d0[4] type', type(d[4]).__name__, 'len', len(d[4]))
    if i in (1, 2):
        print('first item len', len(d[0]), 'repr first 10', repr(d[0][:10])[:400])
        print('types first item', [type(x).__name__ for x in d[0][:15]])
