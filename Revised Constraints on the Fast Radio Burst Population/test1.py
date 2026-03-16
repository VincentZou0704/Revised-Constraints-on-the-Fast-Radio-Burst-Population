import matplotlib.pyplot as plt
import numpy as np

# 创建一个图形
fig, ax = plt.subplots()

# 在范围 [0, 100] 上创建数据点
x_data = np.linspace(0, 100, 100)
y_data = np.sin(x_data)  # 这里只是随便选择了一个函数作为示例数据

# 绘制数据
ax.plot(x_data, y_data, color='blue')

# 自定义刻度位置
major_tick_locs = [10, 30, 50, 70, 90]
minor_tick_locs = np.arange(0, 100, 5)

# 设置主要刻度和次要刻度
ax.xaxis.set_major_locator(plt.FixedLocator(major_tick_locs))
ax.xaxis.set_minor_locator(plt.FixedLocator(minor_tick_locs))

# 添加对数刻度
ax2 = ax.twinx()
ax2.set_yscale('log')
ax2.set_ylim(20, 100)
ax2.yaxis.set_major_locator(plt.MultipleLocator(base=10))
ax2.yaxis.set_minor_locator(plt.MultipleLocator(base=2))
ax2.set_ylabel('Log Scale (20-100)')

# 在刻度上添加锯齿形状
ax.tick_params(which='minor', length=10, direction='in', bottom=True)
ax2.tick_params(which='minor', length=10, direction='in', bottom=True)

# 显示图形
plt.show()
