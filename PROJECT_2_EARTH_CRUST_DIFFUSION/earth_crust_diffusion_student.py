"""
学生模板：地壳热扩散数值模拟
文件：earth_crust_diffusion_student.py
重要：函数名称必须与参考答案一致！
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 设置matplotlib支持中文
rcParams['font.sans-serif'] = ['SimHei']  # 或 'Microsoft YaHei'，根据你的系统字体
rcParams['axes.unicode_minus'] = False    # 正确显示负号
def solve_earth_crust_diffusion():
    """
    实现显式差分法求解地壳热扩散问题
    
    返回:
        tuple: (depth_array, temperature_matrix)
        depth_array: 深度坐标数组 (m)
        temperature_matrix: 温度场矩阵 (°C)
    """
    # 1. 设置物理参数
    A = 10.0  # 平均温度 (°C)
    B = 12.0  # 温度振幅 (°C)
    tau = 365.0  # 周期 (天)
    D = 0.1  # 热扩散率 (m²/day)
    T_bottom = 11.0  # 底部温度 (°C)
    max_depth = 20.0  # 最大深度 (m)
    total_years = 10  # 模拟总年数
    days_per_year = 365
    total_days = total_years * days_per_year
    
    # 2. 设置网格参数 (需要满足稳定性条件)
    dz = 0.1  # 空间步长 (m)
    # 调整时间步长以满足稳定性条件
    dt = 0.004  # 新的时间步长 (天) - 满足 dt <= dz²/(2D) = 0.05
    
    # 计算网格点数
    n_depth = int(max_depth / dz) + 1
    n_time = int(total_days / dt) + 1
    
    # 创建网格
    depth_array = np.linspace(0, max_depth, n_depth)
    time_array = np.linspace(0, total_days, n_time)
    
    # 初始化温度场
    T = np.zeros((n_time, n_depth))
    
    # 3. 设置初始条件 (线性分布)
    T[0, :] = np.linspace(A, T_bottom, n_depth)
    
    # 计算稳定性系数
    r = D * dt / (dz ** 2)
    print(f"稳定性系数 r = {r:.3f}")
    if r > 0.5:
        raise ValueError(f"稳定性条件不满足: r = {r:.3f} > 0.5")
    
    # 4. 时间推进
    for n in range(n_time - 1):
        current_time = time_array[n]
        
        # 应用边界条件
        # 上边界: 时变温度
        T[n+1, 0] = A + B * np.sin(2 * np.pi * current_time / tau)
        # 下边界: 固定温度
        T[n+1, -1] = T_bottom
        
        # 内部点: 显式差分
        for i in range(1, n_depth - 1):
            T[n+1, i] = T[n, i] + r * (T[n, i+1] - 2*T[n, i] + T[n, i-1])
    
    return depth_array, T, tau  # 现在也返回tau

def analyze_results(depth_array, T, tau):
    """分析并可视化结果
    参数:
        depth_array: 深度数组
        T: 温度场矩阵
        tau: 周期 (天)
    """
    # 参数设置
    dt = 0.004  # 与solve函数中一致
    days_per_year = 365
    n_time = T.shape[0]
    time_array = np.linspace(0, 10*365, n_time)
    
    # 由于数据量很大，我们每隔100个点采样一次以加快绘图
    sample_step = 100
    
    # 1. 长期演化分析 - 选择几个深度点观察温度随时间变化
    depths_to_analyze = [0, 1, 5, 10, 15, 20]  # 地表, 1m, 5m, 10m, 15m, 20m
    depth_indices = [np.argmin(np.abs(depth_array - d)) for d in depths_to_analyze]
    
    plt.figure(figsize=(12, 6))
    for i, idx in enumerate(depth_indices):
        plt.plot(time_array[::sample_step]/days_per_year, T[::sample_step, idx], 
                label=f'{depth_array[idx]:.1f}m')
    plt.xlabel('时间 (年)')
    plt.ylabel('温度 (°C)')
    plt.title('不同深度温度随时间变化')
    plt.legend()
    plt.grid()
    plt.show()
    
    # 2. 计算振幅衰减和相位延迟
    # 取最后一年的数据
    last_year_days = int(365 / dt)
    last_year = T[-last_year_days:]
    
    # 计算每个深度的振幅和相位
    amplitudes = []
    phase_shifts = []
    for i in range(len(depth_array)):
        temp_series = last_year[:, i]
        # 使用FFT计算振幅和相位
        fft = np.fft.fft(temp_series - np.mean(temp_series))
        freq = np.fft.fftfreq(len(temp_series), d=dt)
        idx = np.argmax(np.abs(fft)[1:]) + 1  # 忽略直流分量
        
        amplitude = np.abs(fft[idx]) / (len(temp_series)/2)
        phase = np.angle(fft[idx])
        
        amplitudes.append(amplitude)
        phase_shifts.append(phase)
    
    # 振幅衰减
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(depth_array, amplitudes)
    plt.xlabel('深度 (m)')
    plt.ylabel('温度振幅 (°C)')
    plt.title('温度振幅随深度变化')
    plt.grid()
    
    # 相位延迟
    plt.subplot(1, 2, 2)
    plt.plot(depth_array, np.array(phase_shifts) * tau / (2*np.pi))
    plt.xlabel('深度 (m)')
    plt.ylabel('相位延迟 (天)')
    plt.title('温度相位延迟随深度变化')
    plt.grid()
    
    plt.tight_layout()
    plt.show()
    
    # 3. 季节性温度轮廓可视化 (第10年的四季)
    seasons = {
        '春季': 365*9 + 80,  # 第10年3月1日左右
        '夏季': 365*9 + 172, # 第10年6月21日左右
        '秋季': 365*9 + 264, # 第10年9月21日左右
        '冬季': 365*9 + 355  # 第10年12月21日左右
    }
    
    plt.figure(figsize=(10, 6))
    for season, day in seasons.items():
        time_idx = int(day / dt)  # 转换为时间索引
        plt.plot(T[time_idx, :], depth_array, label=season)
    
    plt.gca().invert_yaxis()
    plt.xlabel('温度 (°C)')
    plt.ylabel('深度 (m)')
    plt.title('第10年四季温度随深度分布')
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # 运行模拟
    depth, T, tau = solve_earth_crust_diffusion()  # 现在接收tau
    print(f"计算完成，温度场形状: {T.shape}")
    
    # 分析结果
    analyze_results(depth, T, tau)  # 传递tau到分析函数
