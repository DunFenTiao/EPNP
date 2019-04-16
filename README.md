# Read Me

@last update: 4.12

## HomeWork 4.13-4.20

1. EPNP
According to the paper:[EPNP](https://icwww.epfl.ch/~lepetit/papers/lepetit_ijcv08.pdf),complete the EPNP functions. 
> Lepetit V, Moreno-Noguer F, Fua P. Epnp: An accurate o (n) solution to the pnp problem[J]. International journal of computer vision, 2009, 81(2): 155.

2. Including

- Compiling with Cmake tool chain.
- Available Libray: Eigen , openCV
- Unit Testing: google-test

3. to be Submitted

- Project EPNP
- summary report(PPT) 

✅ 环境配置 cmake,opencv343,eigen
✅ feature extraction, matching
✅ opencv solver,eigen solver
todo: gtest, PPT summary
![img](https://wx3.sinaimg.cn/mw690/c7716318ly1g2172ibxq1j20ze0dx7wh.jpg)


> result: 从打印结果R，T看，两种解法结果基本相同
-- Max dist : 95.000000 
-- Min dist : 7.000000 
一共找到了81组匹配点
3d-2d pairs: 77

---using opencv solvePNP---

R1: [0.9978654469053291, -0.05171629787598182, 0.039874483149384;
 0.0505533717792531, 0.9982813811698801, 0.02964187260119196;
 -0.04133892202484731, -0.02756281087914287, 0.9987649297919227]
t1:[-0.1266609334454444;
 -0.0111141371771197;
 0.05673412814657741]

 ---using egien to solve EPNP---

R1:   0.997854 -0.0518631  0.0399726
 0.0506961   0.998273  0.0296762
-0.0414427 -0.0275861    0.99876
t1: -0.127002
-0.0109942
 0.0567606

# EPNP 原理总结

## 选取控制点

世界坐标系下3D参考点
$$p_i, i=1,...,n.$$
四个控制点
$$c_j, j=1,..,4.$$
每个参考点都可以用控制点来表示,$a_{ij}$为均一重心坐标系（homogeneous barycentric coordinates）
$$p_i^w= \sum^4_{j=1}a_{ij}c_i^w,with\sum^4_{j=1}a_{ij}=1$$
相机坐标系下的3D参考点同理也可以这样表示
$$p^c_i=\sum^4_{j=1}a_{ij}c^c_j.$$
选择参考点们的重心作为第一个控制点，剩下三个选择和主轴方向一致的点，可以增加解的稳定程度。
![此处输入图片的描述][1]
## 计算权重和的特征向量
3D参考点$\bf{p}_i$在像素平面的投影为$\bf{u_i}$,$\bf{A}$为相机内参矩阵，$w_i$是尺度系数
$$
w_i\begin{bmatrix}\bf{u}_i\\1
\end{bmatrix} 
= Ap_i^c 
= A\sum^4_{j=1}a_{ij}c^c_j.
$$
上式中，3D控制点用$\begin{bmatrix}x_j^c\\y_j^c\\z_j^c
\end{bmatrix} $来表示，控制点在像素平面上投影用$\begin{bmatrix}u_i\\v_i\\1
\end{bmatrix} $来表示，$A$用$\begin{bmatrix}f_u & 0 & u_c\\0 & f_v & v_c\\0 & 0 & 1
\end{bmatrix} $来表示，则有
$$
w_i\begin{bmatrix}u_i\\v_i\\1
\end{bmatrix} 
= \begin{bmatrix}f_u & 0 & u_c\\0 & f_v & v_c\\0 & 0 & 1
\end{bmatrix} 
\sum^4_{j=1}a_{ij}
\begin{bmatrix}x_j^c\\y_j^c\\z_j^c
\end{bmatrix} 
.
$$
上式可以消去$w_i$,化为两个线性的方程组
$$
\sum^4_{j=1}a_{ij}f_ux^c_j+a_{ij}(u_c-u_i)z^c_j=0,\\
\sum^4_{j=1}a_{ij}f_vy^c_j+a_{ij}(v_c-v_i)z^c_j=0,
$$
这两个方程组只关于$\bf{x}=[c_1^{c^T},c_2^{c^T},c_3^{c^T},c_4^{c^T}]^t$,
可以统一为
$$
\bf{Mx=0}\\
\bf{x}=\sum^N_{j=1}\beta_i\bf{v}_i
$$
x属于M的右零空间，令$\bf{v}_i$为矩阵M的右奇异向量，$\bf{v}_i$可以通过求解$MTM$的零空间特征值得到。

## 选择合适的线性组合
特征值的个数不是固定的，EPNP计算四种情况1-4，选择重投影误差最小的那种情况，其中dist是像素平面上点，和重投影点齐次坐标的2D距离。
$$
res=\sum_idist^2(\bf{A[R|t]}
\begin{bmatrix}\bf{p}_i^w\\1
\end{bmatrix} ，\bf{u}_i)
$$
在相机坐标系，和世界坐标系下，控制点之间的距离应该是相同的
$$
\begin{Vmatrix}c^c_i-c^c_j\end{Vmatrix}^2=\begin{Vmatrix}c^w_i-c^w_j\end{Vmatrix}^2$$

$$c^c_i=\sum^4_{j=1}\beta_jv^{[i]}_j.
$$
四个控制点，可以有6组这样的方程，

- 情况一： N=1，未知数个数1
$$ \bf{x} =\beta{\bf{v}}$$
- 情况二： N=2，未知数个数3
 $$ \bf{x} =\beta_1{\bf{v_1}}+\beta_2{\bf{v_2}}$$
- 情况三： N=3，未知数个数6
$$\bf{L}\beta=\bf{\rho}\\
\bf{L}=[\bf{v_1},\bf{v_2},\bf{v_3}]\\
\beta=[\beta_{11},\beta_{12},\beta_{13},\beta_{22},\beta_{23},\beta_{33}]$$
- 情况四： N=4，未知数个数10


##高斯牛顿优化

$$
Error(\bf{\beta})=\sum_{(i,j)s.t. i<j}(\begin{Vmatrix}c^c_i-c^c_j\end{Vmatrix}^2-\begin{Vmatrix}c^w_i-c^w_j\end{Vmatrix}^2)\\
c^c_i=\sum^4_{j=1}\beta_jv^{[i]}_j.
$$

## 求相机姿态

1. 控制点在相机参考系下坐标
$$c^c_i=\sum^4_{j=1}\beta_jv^{[i]}_j$$
2. 参考点在相机坐标系下坐标
$$p^c_i=\sum^4_{j=1}a_{ij}c^c_j.$$
3. 计算世界坐标系下参考点重心
$$
p_0^w=\frac1n\sum^n_{i=1}p_i^w
$$
4. 计算相机坐标系下参考点重心
$$
p_0^c=\frac1n\sum^n_{i=1}p_i^c
$$
5. 计算H
$$
H= B^TA
$$
6. 对于H的SVD分解
$$
H= UΣV^T
$$
7. 计算位姿中旋转
$$
R= UV^T
$$
7. 计算位姿中平移
$$
t= p^c_0-Rp_0^w
$$

  [1]: https://wx3.sinaimg.cn/mw690/c7716318ly1g24kwqwl4vj20i00ao3zr.jpg

