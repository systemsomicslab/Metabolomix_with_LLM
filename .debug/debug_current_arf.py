import os
import io
import msgpack
import lz4.block
from metabolomix.arf.reader import deserialize_lz4_packed_msgpack

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
path = r'data\AlignmentResult_2026_05_14_15_16_23_PeakProperties.arf'
with open(path, 'rb') as f:
    data = f.read()
all_data = deserialize_lz4_packed_msgpack(data)
print('top count', len(all_data))
for i, d in enumerate(all_data):
    print(f'--- top {i}, type {type(d).__name__}, len {len(d) if hasattr(d, "__len__") else "N/A"}')
    if i == 0:
        print('first[0]:', d[0])
        print('first[1]:', d[1])
        print('len(first[2:]):', len(d[2:]))
        # Check first few items in first[2:]
        for j in range(min(5, len(d[2:]))):
            item = d[2 + j]
            print(f'  item {j}: type {type(item).__name__}, len {len(item) if hasattr(item, "__len__") else "N/A"}')
            if isinstance(item, list):
                print(f'    item[{j}] content: {item}')
                if len(item) >= 3:
                    print(f'    subitems types: {[type(sub).__name__ for sub in item[:3]]}')