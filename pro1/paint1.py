import matplotlib.pyplot as plt

# 设置全局字体
plt.rcParams.update({'font.family': 'Microsoft YaHei'})  # 这里使用微软雅黑字体，根据需要可以修改成其他字体

# 数据
labels = ['2^8', '2^9', '2^10', '2^11', '2^12']
trivial_algorithm = [0.1325, 0.5490, 2.5479, 15.9271, 87.3309]
cache_optimization = [0.1099, 0.4497, 1.8294, 7.4483, 29.2393]

# 绘制折线图
plt.plot(labels, trivial_algorithm, marker='o', label='平凡算法')
plt.plot(labels, cache_optimization, marker='s', label='逐行访问Cache优化')

# 添加标题和标签
plt.title('不同优化策略的性能对比')
plt.xlabel('数据规模')
plt.ylabel('平均时间(MS)')
plt.legend()

# 显示网格
plt.grid(True)


plt.savefig('martix1.pdf')

# 显示图表
plt.show()