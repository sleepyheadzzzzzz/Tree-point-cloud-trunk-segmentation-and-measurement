import subprocess
import os

#
cc_path = "C:\Program Files\CloudCompare\CloudCompare.exe"

# 
bin_file = "C:\Users\Admin\Desktop\mobile\3_1.bin"

# 
output_dir = "C:\Users\Admin\Desktop\mobile\3_1"

# 
number_of_clouds = 1

# 
for i in range(number_of_clouds):
    output_file = f"{i}.ply"
    # 构建CloudCompare命令
    command = f'"{cc_path}" -SILENT -O "{bin_file}" -C_EXPORT_FMT PLY -PLY_EXPORT_FMT BINARY_BE -NO_TIMESTAMP -SAVE_CLOUDS "{output_file}"'
    # 运行命令
    subprocess.run(command)

