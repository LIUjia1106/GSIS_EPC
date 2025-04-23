import numpy as np
import os

# 数据参数
num_kpoints = 512
num_e_bands = 3
num_qpoints = 1728
num_p_bands = 6

# 初始化矩阵
C_matrix = np.zeros((num_kpoints * num_e_bands, num_qpoints * num_p_bands))

# 文件夹路径
folder_path = './input/ep_matrix/'  # 替换为实际文件夹路径

# 遍历所有文件
for k in range(1, num_kpoints + 1):
    file_name = f'ep_matrix_{k}_200.000007464979.txt'  # 文件命名规则
    file_path = os.path.join(folder_path, file_name)
    
    if not os.path.exists(file_path):
        print(f"File {file_name} not found. Skipping.")
        continue
    
    # 读取文件
    data = np.loadtxt(file_path)
    
    # 解析每行数据并存储到矩阵
    for row in data:
        qpoint_index = int(row[3]) - 1  # 声子 qpoint index (0-based)
        phonon_band_index = int(row[4]) - 1  # 声子 band index (0-based)
        kpoint_index = int(row[6]) - 1  # 电子 kpoint index (0-based)
        electron_band_index = int(row[7]) - 1  # 电子 band index (0-based)
        C_value = row[-1]  # C 值
        
        # 矩阵索引
        row_idx = kpoint_index * num_e_bands + electron_band_index
        col_idx = qpoint_index * num_p_bands + phonon_band_index
        
        # 存储 C 值
        C_matrix[row_idx, col_idx] = C_value

# 保存矩阵为ASCI文件
np.savetxt('C_matrix.txt', C_matrix, fmt='%.6e')
# 保存为 NumPy 文件，便于 Fortran 调用
# np.save('C_matrix.npy', C_matrix)  

# 直接以 Fortran 风格的二进制文件保存
# C_matrix.flatten(order='F').tofile('C_matrix.dat')
# 在Fortran中的调用方式
# real(8), dimension(512*3, 1728*6) :: C_matrix
# open(unit=10, file='C_matrix.dat', form='unformatted', access='stream')
# read(10) C_matrix
# close(10)


