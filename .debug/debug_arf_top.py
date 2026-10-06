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
print('top_count', len(all_data))
for i,d in enumerate(all_data):
    print('--- d', i, type(d).__name__, 'len', len(d))
    if isinstance(d, list) and len(d) < 14:
        for j,item in enumerate(d):
            print('  item', j, type(item).__name__, 'len', len(item) if hasattr(item,'__len__') else None)
            if isinstance(item, list) and len(item) > 14:
                print('   first 5 types', [type(v).__name__ for v in item[:5]])
                print('   repr first 200', repr(item[:3])[:200])
    elif isinstance(d, list) and len(d) >= 14:
        print('  first 5 types', [type(v).__name__ for v in d[:5]])
        print('  len d[2]', len(d[2]) if isinstance(d[2], list) else 'NA')
        print('  len d[3]', len(d[3]) if isinstance(d[3], list) else 'NA')
        print('  len d[4]', len(d[4]) if isinstance(d[4], list) else 'NA')
        print('  len d[10]', len(d[10]) if len(d)>10 and isinstance(d[10], list) else 'NA')
        print('  len d[16]', len(d[16]) if len(d)>16 and isinstance(d[16], list) else 'NA')
