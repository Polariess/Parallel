import matplotlib.pyplot as plt

# 设置全局字体
plt.rcParams.update({'font.family': 'Microsoft YaHei'})  # 这里使用微软雅黑字体，根据需要可以修改成其他字体

# 数据
labels = ['2^8', '2^9', '2^10', '2^11', '2^12']
speedup_ratio = [1.21, 1.22, 1.39, 2.14, 2.99]

# 绘制折线图
plt.plot(labels, speedup_ratio, linestyle='--',marker='o')

# 添加标题和标签
plt.title('加速比增长曲线')
plt.xlabel('数据规模')
plt.ylabel('加速比')
plt.grid(True)

plt.savefig('martix2.pdf')

# 显示图表
plt.show()
