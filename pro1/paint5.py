import matplotlib.pyplot as plt

# 设置全局字体为微软雅黑
plt.rcParams.update({'font.family': 'Microsoft YaHei'})

# 数据
labels = ['$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$']
trivial_algorithm = [1.4728, 1.4409, 1.4493, 1.4368, 1.4409]
double_linked = [2.6738, 2.4876, 2.5575, 2.5907, 2.5189]
recursive_call = [5.1813, 5.0251, 5.1282, 5.1546, 5.1020]

# 绘制折线图
plt.plot(labels, trivial_algorithm, marker='o', label='平凡算法')
plt.plot(labels, double_linked, marker='s', label='双链路式')
plt.plot(labels, recursive_call, marker='^', label='类递归调用')

# 添加标题和标签
plt.title('不同算法与数据规模下的IPC指标')
plt.xlabel('数据规模')
plt.ylabel('IPC指标')
plt.legend()

# 显示网格
plt.grid(True)

# 保存为pdf文件
plt.savefig('ipc_indicator.pdf')

# 显示图表
plt.show()
