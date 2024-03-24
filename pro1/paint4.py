import matplotlib.pyplot as plt

# 设置全局字体
plt.rcParams.update({'font.family': 'Microsoft YaHei'})  # 设置为微软雅黑字体，根据需要可以修改成其他字体

# 数据
labels = ['$2^{8}$', '$2^{9}$', '$2^{10}$', '$2^{11}$', '$2^{12}$']
trivial_algorithm = [96.91, 89.58, 66.07, 66.10, 66.42]
cache_optimized = [98.51, 98.37, 98.78, 99.11, 98.53]

# 绘制折线图
plt.plot(labels, trivial_algorithm, marker='o', label='平凡算法')
plt.plot(labels, cache_optimized, marker='s', label='逐行访问Cache优化')

# 添加标题和标签
plt.title('不同算法与数据规模下的Cache命中率')
plt.xlabel('数据规模')
plt.ylabel('Cache命中率 (%)')
plt.legend()

# 显示网格
plt.grid(True)

# 设置y轴最低为0
plt.ylim(0)

# 保存为pdf文件
plt.savefig('cache_hit_rate.pdf')

# 显示图表
plt.show()
