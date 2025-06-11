"""学生模板：量子隧穿效应
文件：quantum_tunneling_student.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

class QuantumTunnelingSolver:
    """量子隧穿求解器类
    
    该类实现了一维含时薛定谔方程的数值求解，用于模拟量子粒子的隧穿效应。
    使用变形的Crank-Nicolson方法进行时间演化，确保数值稳定性和概率守恒。
    """
    
    def __init__(self, Nx=220, Nt=300, x0=40, k0=0.5, d=10, barrier_width=3, barrier_height=1.0):
        """初始化量子隧穿求解器
        
        参数:
            Nx (int): 空间网格点数，默认220
            Nt (int): 时间步数，默认300
            x0 (float): 初始波包中心位置，默认40
            k0 (float): 初始波包动量(波数)，默认0.5
            d (float): 初始波包宽度参数，默认10
            barrier_width (int): 势垒宽度，默认3
            barrier_height (float): 势垒高度，默认1.0
        """
        # TODO: 初始化所有参数
        # 提示：需要设置空间网格、势垒参数，并初始化波函数和系数矩阵
        self.Nx = Nx
        self.Nt = Nt
        self.x0 = x0
        self.k0 = k0
        self.d = d
        self.barrier_width = int(barrier_width)  # 确保是整数
        self.barrier_height = barrier_height
        
        # TODO: 创建空间网格
        self.x = np.arange(self.Nx)  # 应该是 np.arange(self.Nx)
        
        # TODO: 设置势垒
        self.V = self.setup_potential()  # 调用 setup_potential() 方法

        # TODO: 初始化波函数矩阵和系数矩阵
        self.C = np.zeros((self.Nx, self.Nt), dtype=np.complex128)  # 复数矩阵，形状为 (Nx, Nt)
        self.B = np.zeros((self.Nx, self.Nt), dtype=np.complex128)  # 复数矩阵，形状为 (Nx, Nt)

    def wavefun(self, x):
        """高斯波包函数
        
        参数:
            x (np.ndarray): 空间坐标数组
            
        返回:
            np.ndarray: 初始波函数值
            
        数学公式:
            ψ(x,0) = exp(ik₀x) * exp(-(x-x₀)²ln10(2)/d²)
            
        物理意义:
            描述一个在x₀位置、具有动量k₀、宽度为d的高斯波包
        """
        # TODO: 实现高斯波包函数
        # 提示：包含动量项 exp(ik₀x) 和高斯包络 exp(-(x-x₀)²ln10(2)/d²)
        return np.exp(self.k0 * 1j * x) * np.exp(-((x - self.x0) ** 2) * np.log10(2) / (self.d ** 2))

    def setup_potential(self):
        """设置势垒函数
        
        返回:
            np.ndarray: 势垒数组
            
        说明:
            在空间网格中间位置创建矩形势垒
            势垒位置：从 Nx//2 到 Nx//2+barrier_width
            势垒高度：barrier_height
        """
        # TODO: 创建势垒数组
        # 提示：
        # 1. 初始化全零数组
        # 2. 在中间位置设置势垒高度
        # 3. 注意barrier_width必须是整数
        V = np.zeros(self.Nx)
        V[self.Nx // 2:self.Nx // 2 + self.barrier_width] = self.barrier_height
        return V
    
    def build_coefficient_matrix(self):
        """构建变形的Crank-Nicolson格式的系数矩阵
        
        返回:
            np.ndarray: 系数矩阵A
            
        数学原理:
            对于dt=1, dx=1的情况，哈密顿矩阵的对角元素为: -2+2j-V
            非对角元素为1（表示动能项的有限差分）
            
        矩阵结构:
            三对角矩阵，主对角线为 -2+2j-V[i]，上下对角线为1
        """
        # TODO: 构建系数矩阵
        # 提示：
        # 1. 使用 np.diag() 创建三对角矩阵
        # 2. 主对角线：-2+2j-self.V
        # 3. 上对角线和下对角线：全1数组
        A = np.diag(-2 + 2j - self.V) + np.diag(np.ones(self.Nx - 1), 1) + np.diag(np.ones(self.Nx - 1), -1)
        return A

    def solve_schrodinger(self):
        """求解一维含时薛定谔方程
        
        使用Crank-Nicolson方法进行时间演化
        
        返回:
            tuple: (x, V, B, C) - 空间网格, 势垒, 波函数矩阵, chi矩阵
            
        数值方法:
            Crank-Nicolson隐式格式，具有二阶精度和无条件稳定性
            时间演化公式：C[:,t+1] = 4j * solve(A, B[:,t])
                         B[:,t+1] = C[:,t+1] - B[:,t]
        """
        # TODO: 实现薛定谔方程求解
        # 提示：
        # 1. 构建系数矩阵A
        # 2. 设置初始波函数 B[:,0] = wavefun(x)
        # 3. 对初始波函数进行归一化
        # 4. 时间循环：使用线性方程组求解进行时间演化
        A = self.build_coefficient_matrix()  # 构建系数矩阵
        self.B[:, 0] = self.wavefun(self.x)  # 设置初始波函数
        self.B[:, 0] /= np.sqrt(np.sum(np.abs(self.B[:, 0])**2))  # 归一化初始波函数
        for t in range(self.Nt - 1):
            self.C[:, t + 1] = 4j * np.linalg.solve(A, self.B[:, t])
            self.B[:, t + 1] = self.C[:, t + 1] - self.B[:, t]
        return self.x, self.V, self.B, self.C

    def calculate_coefficients(self):
        """计算透射和反射系数
        
        返回:
            tuple: (T, R) - 透射系数和反射系数
            
        物理意义:
            透射系数T：粒子穿过势垒的概率
            反射系数R：粒子被势垒反射的概率
            应满足：T + R ≈ 1（概率守恒）
            
        计算方法:
            T = ∫|ψ(x>barrier)|²dx / ∫|ψ(x)|²dx
            R = ∫|ψ(x<barrier)|²dx / ∫|ψ(x)|²dx
        """
        # TODO: 计算透射和反射系数
        # 提示：
        # 1. 确定势垒位置
        # 2. 计算透射区域的概率（势垒右侧）
        # 3. 计算反射区域的概率（势垒左侧）
        # 4. 归一化处理
        barrier_position = self.Nx // 2
        transmitted_prob = np.sum(np.abs(self.B[barrier_position + self.barrier_width:, -1])**2)
        reflected_prob = np.sum(np.abs(self.B[:barrier_position, -1])**2)
        total_prob = np.sum(np.abs(self.B[:, -1])**2)
        return transmitted_prob / total_prob, reflected_prob / total_prob

    def plot_evolution(self, time_indices=None):
        """绘制波函数演化图
        
        参数:
            time_indices (list): 要绘制的时间索引列表，默认为[0, Nt//4, Nt//2, 3*Nt//4, Nt-1]
            
        功能:
            在多个子图中显示不同时刻的波函数概率密度和势垒
        """
        # TODO: 实现波函数演化绘图
        # 提示：
        # 1. 设置默认时间索引
        # 2. 创建子图布局
        # 3. 绘制概率密度 |ψ|²
        # 4. 绘制势垒
        # 5. 添加标题和标签
        if time_indices is None:
            Nt = self.B.shape[1]
            time_indices = [0, Nt//4, Nt//2, 3*Nt//4, Nt-1]
        
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        axes = axes.flatten()
        
        # Add overall title with barrier parameters
        fig.suptitle(f'Quantum Tunneling Evolution - Barrier Width: {self.barrier_width}, Barrier Height: {self.barrier_height}', 
                     fontsize=14, fontweight='bold')
        
        for i, t_idx in enumerate(time_indices):
            if i < len(axes):
                ax = axes[i]
                
                # Plot probability density
                prob_density = np.abs(self.B[:, t_idx])**2
                ax.plot(self.x, prob_density, 'b-', linewidth=2, 
                       label=f'|ψ|² at t={t_idx}')
                
                # Plot potential
                ax.plot(self.x, self.V, 'k-', linewidth=2, 
                       label=f'Barrier (Width={self.barrier_width}, Height={self.barrier_height})')
                
                ax.set_xlabel('Position')
                ax.set_ylabel('Probability Density')
                ax.set_title(f'Time step: {t_idx}')
                ax.legend()
                ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    def create_animation(self, interval=20):
        """创建波包演化动画
        
        参数:
            interval (int): 动画帧间隔(毫秒)，默认20
            
        返回:
            matplotlib.animation.FuncAnimation: 动画对象
            
        功能:
            实时显示波包在势垒附近的演化过程
        """
        # TODO: 创建动画
        # 提示：
        # 1. 设置图形和坐标轴
        # 2. 创建线条对象
        # 3. 定义动画更新函数
        # 4. 使用 FuncAnimation 创建动画
        Nx, Nt = self.B.shape
        fig = plt.figure(figsize=(10, 6))
        plt.axis([0, Nx, 0, np.max(self.V)*1.1])
        
        # Add title with barrier parameters
        plt.title(f'Quantum Tunneling Animation - Barrier Width: {self.barrier_width}, Barrier Height: {self.barrier_height}', 
                 fontsize=12, fontweight='bold')
        plt.xlabel('Position')
        plt.ylabel('Probability Density / Potential')
        
        myline, = plt.plot([], [], 'r', lw=2, label='|ψ|²')
        myline1, = plt.plot(self.x, self.V, 'k', lw=2, 
                           label=f'Barrier (Width={self.barrier_width}, Height={self.barrier_height})')
        
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        def animate(i):
            myline.set_data(self.x, np.abs(self.B[:, i]))
            myline1.set_data(self.x, self.V)
            return myline, myline1
        
        anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=interval)
        return anim

    def verify_probability_conservation(self):
        """验证概率守恒

        返回:
            np.ndarray: 每个时间步的总概率

        物理原理:
            量子力学中概率必须守恒：∫|ψ(x,t)|²dx = 常数
            数值计算中应该保持在1附近
        """
        # TODO: 验证概率守恒
        # 提示：
        # 1. 计算每个时间步的总概率
        # 2. 考虑空间步长dx的影响
        # 3. 返回概率数组用于分析
        total_prob = np.zeros(self.Nt)
        for t in range(self.Nt):
            total_prob[t] = np.sum(np.abs(self.B[:, t])**2)
        return total_prob

    def demonstrate(self):
        """演示量子隧穿效应
        
        功能:
            1. 求解薛定谔方程
            2. 计算并显示透射和反射系数
            3. 绘制波函数演化图
            4. 验证概率守恒
            5. 创建并显示动画
            
        返回:
            animation对象
        """
        # TODO: 实现完整的演示流程
        # 提示：
        # 1. 打印开始信息
        # 2. 调用solve_schrodinger()
        # 3. 计算并显示系数
        # 4. 绘制演化图
        # 5. 验证概率守恒
        # 6. 创建动画
        self.solve_schrodinger()
        T, R = self.calculate_coefficients()
        print(f"\n势垒宽度:{self.barrier_width}, 势垒高度:{self.barrier_height} 结果")
        print(f"Transmission coefficient: {T:.4f}")
        print(f"Reflection coefficient: {R:.4f}")
        print(f"Total (T + R): {T + R:.4f}")
        # Plot evolution
        print("\nPlotting wave function evolution...")
        self.plot_evolution()
        
        # Check probability conservation
        total_prob = self.verify_probability_conservation()
        print(f"\nProbability conservation:")
        print(f"Initial probability: {total_prob[0]:.6f}")
        print(f"Final probability: {total_prob[-1]:.6f}")
        print(f"Relative change: {abs(total_prob[-1] - total_prob[0])/total_prob[0]*100:.4f}%")
        
        # Create animation
        print("\nCreating animation...")
        anim = self.create_animation()
        plt.show()
        
        return anim


def demonstrate_quantum_tunneling():
    """便捷的演示函数
    
    创建默认参数的求解器并运行演示
    
    返回:
        animation对象
    """
    # TODO: 创建求解器实例并调用demonstrate方法
    solver = QuantumTunnelingSolver()
    return solver.demonstrate()
    


if __name__ == "__main__":
    # 运行演示
    barrier_width = 3
    barrier_height = 1.0
    solver = QuantumTunnelingSolver(barrier_width=barrier_width, barrier_height=barrier_height)
    animation = solver.demonstrate()
    
