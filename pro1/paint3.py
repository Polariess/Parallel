import matplotlib.pyplot as plt

# 设置全局字体
plt.rcParams.update({'font.family': 'Microsoft YaHei'})  # 这里使用微软雅黑字体，根据需要可以修改成其他字体

# 数据
labels = ['$2^{16}$', '$2^{17}$', '$2^{18}$', '$2^{19}$', '$2^{20}$']
trivial_algorithm = [0.1076, 0.2486, 0.4952, 0.9727, 1.9903]
double_linked = [0.0574, 0.1306, 0.2593, 0.5042, 1.0380]
recursive_call = [0.0823, 0.1792, 0.3514, 0.7038, 1.4001]

# 绘制折线图
plt.plot(labels, trivial_algorithm, marker='o', label='平凡算法')
plt.plot(labels, double_linked, marker='s', label='双链路式')
plt.plot(labels, recursive_call, marker='D', label='类递归调用')

# 添加标题和标签
plt.title('不同算法与数据规模下的性能对比')
plt.xlabel('数据规模')
plt.ylabel('平均时间(MS)')
plt.legend()

# 显示网格
plt.grid(True)

# 保存为pdf文件
plt.savefig('add1.pdf')

# 显示图表
plt.show()
