
"""
学生模板：铝棒热传导问题
文件：heat_diffusion_student.py
重要：函数名称必须与参考答案一致！
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 物理参数
K = 237       # 热导率 (W/m/K)
C = 900       # 比热容 (J/kg/K)
rho = 2700    # 密度 (kg/m^3)
D = K/(C*rho) # 热扩散系数
L = 1         # 铝棒长度 (m)
dx = 0.01     # 空间步长 (m)
dt = 0.5      # 时间步长 (s)
Nx = int(L/dx) + 1 # 空间格点数
Nt = 2000     # 时间步数

def basic_heat_diffusion():
    """
    任务1: 基本热传导模拟
    
    返回:
        np.ndarray: 温度分布数组
    """
    r = D * dt / dx**2# 稳定性参数
    print(f"[任务1] 稳定性参数 r = {r:.4f}")
    u = np.zeros((Nx, Nt))# 创建温度数组，初始全为0
    u[1:-1, 0] = 100.0  # 设置初始条件，内部点温度为100K
    for n in range(Nt-1):# 时间步循环
        # 显式有限差分法更新内部点温度
        u[1:-1, n+1] = (1-2*r)*u[1:-1, n] + r*(u[2:, n] + u[:-2, n])
        u[0, n+1] = 0.0# 左端边界条件
        u[-1, n+1] = 0.0# 右端边界条件
    return u # 返回温度分布数组
def analytical_solution(n_terms=100):
    """
    任务2: 解析解函数
    
    参数:
        n_terms (int): 傅里叶级数项数
    
    返回:
        np.ndarray: 解析解温度分布
    """
    x = np.linspace(0, L, Nx) # 生成空间坐标数组，Nx个点，范围[0, L]
    t = np.arange(Nt) * dt# 生成时间数组，Nt个点，步长为dt
    T0 = 100 # 初始温度T0=100K
    u = np.zeros((Nx, Nt))# 创建温度分布数组，初始全为0
    for n in range(n_terms):# 对前n_terms项傅里叶级数求和
        k_n = (2*n+1)*np.pi/L# 计算第n项的波数k_n
        coef = 4*T0/((2*n+1)*np.pi) # 计算第n项的系数
        sin_knx = np.sin(k_n * x[:, None])# 计算sin(k_n x)，对所有x，形状(Nx,1)
        exp_term = np.exp(-k_n**2 * D * t[None, :])# 计算指数衰减项，对所有t，形状(1,Nt)
        u += coef * sin_knx * exp_term# 叠加第n项的贡献到总温度分布
    return u# 返回解析解温度分布数组

def stability_analysis():
    """
    任务3: 数值解稳定性分析
    """
    dx_ = 0.01# 空间步长
    dt_ = 0.6  # 时间步长使r>0.5
    r = D * dt_ / dx_**2#稳定性参数
    print(f"[任务3] 尝试不稳定步长 dt={dt_}, r={r:.4f} (应>0.5)")
    Nx_ = int(L/dx_) + 1# 计算空间格点数
    Nt_ = 2000# 设置时间步数
    u = np.zeros((Nx_, Nt_))# 创建温度数组，初始全为0
    u[1:-1, 0] = 100.0# 设置初始条件，内部点温度为100K
    for n in range(Nt_-1):
        u[1:-1, n+1] = (1-2*r)*u[1:-1, n] + r*(u[2:, n] + u[:-2, n])
        u[0, n+1] = 0.0
        u[-1, n+1] = 0.0
    print(f"[任务3] 数值解最大值: {np.max(u)}")
    print(f"[任务3] 数值解最小值: {np.min(u)}")
    print("[任务3] 若出现极大/极小值或震荡，说明已不稳定。")
    # 可视化（可选）
    plot_3d_solution(u, dx_, dt_, Nt_, "任务3：不稳定步长下的温度分布")

def different_initial_condition():
    """
    任务4: 不同初始条件模拟
    
    返回:
        np.ndarray: 温度分布数组
    """
    Nt_diff = 1000# 设置时间步数为1000
    u = np.zeros((Nx, Nt_diff))# 创建温度数组，初始全为0
    x = np.linspace(0, L, Nx)# 生成空间坐标数组
    u[:, 0] = np.where(x < 0.5, 100.0, 50.0)# 左半段100K，右半段50K
    u[0, 0] = 0.0
    u[-1, 0] = 0.0
    r = D * dt / dx**2
    for n in range(Nt_diff-1):
        u[1:-1, n+1] = (1-2*r)*u[1:-1, n] + r*(u[2:, n] + u[:-2, n])
        u[0, n+1] = 0.0
        u[-1, n+1] = 0.0
    return u# 返回温度分布数组

def heat_diffusion_with_cooling():
    """
    任务5: 包含牛顿冷却定律的热传导
    """
    h = 0.1  # 冷却系数 (s^-1)
    u = np.zeros((Nx, Nt))# 创建温度数组，初始全为0
    u[1:-1, 0] = 100.0# 设置初始条件，内部点温度为100K
    r = D * dt / dx**2# 计算稳定性参数 r
    for n in range(Nt-1):# 时间步循环
        u[1:-1, n+1] = (1-2*r-h*dt)*u[1:-1, n] + r*(u[2:, n] + u[:-2, n])
        u[0, n+1] = 0.0
        u[-1, n+1] = 0.0
    # 按测试要求，不返回任何内容
    return None

def plot_3d_solution(u, dx, dt, Nt, title):
    """
    绘制3D温度分布图
    
    参数:
        u (np.ndarray): 温度分布数组
        dx (float): 空间步长
        dt (float): 时间步长
        Nt (int): 时间步数
        title (str): 图表标题
    
    返回:
        None
    
    示例:
        >>> u = np.zeros((100, 200))
        >>> plot_3d_solution(u, 0.01, 0.5, 200, "示例")
    """
    x = np.linspace(0, dx*(u.shape[0]-1), u.shape[0]) # 生成空间坐标数组
    t = np.arange(Nt) * dt
    X, T = np.meshgrid(x, t, indexing='ij')
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, u, cmap='viridis')
    ax.set_xlabel('position x (m)')
    ax.set_ylabel('time t (s)')
    ax.set_zlabel('temperature T (K)')
    ax.set_title(title)
    fig.colorbar(surf, shrink=0.5, aspect=10, label='temperature (K)')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    """
    主函数 - 演示和测试各任务功能
    
    执行顺序:
    1. 基本热传导模拟
    2. 解析解计算
    3. 数值解稳定性分析
    4. 不同初始条件模拟
    5. 包含冷却效应的热传导
    
    注意:
        学生需要先实现各任务函数才能正常运行
    """
    print("=== 铝棒热传导问题学生实现 ===")
    # 任务1：基本热传导模拟
    u_basic = basic_heat_diffusion()
    plot_3d_solution(u_basic, dx, dt, Nt, "Task 1: Numerical Solution of Basic Heat Diffusion")

    # 任务2：解析解
    u_analytical = analytical_solution(n_terms=100)
    plot_3d_solution(u_analytical, dx, dt, Nt, "Task 2: Analytical Solution")

    # 任务3：数值解稳定性分析
    stability_analysis()

    # 任务4：不同初始条件
    u_diff_init = different_initial_condition()
    plot_3d_solution(u_diff_init, dx, dt, 1000, "Task 4: Different Initial Conditions")

    # 任务5：牛顿冷却
    heat_diffusion_with_cooling()
