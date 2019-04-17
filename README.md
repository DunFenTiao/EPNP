# EPNP README

标签（空格分隔）： SLAM

---

@last update: 4.16

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

✅ 原理总结，画图

❓todo: gtest


![img](https://wx3.sinaimg.cn/mw690/c7716318ly1g2172ibxq1j20ze0dx7wh.jpg)


result: 从打印结果R，T看，两种解法求出的相机姿态基本相同
![r,t][1]
 
 
# EPNP总结
 
## 什么是PNP
根据3D参考点和对应2D的点，估计相机姿态
(the estimation of the pose of a calibrated camera from n 3D-to-2D point correspondences )

- 输入：
N个世界坐标系中的3D参考点坐标
3D点对应的2d参考点坐标
相机内参K

- 输出：
相机位姿
# 原理
用四个控制点的权重和来表示参考点
（express the n 3D points as a weighted sum of four virtual control points 
）

![步骤][3]
![控制点][2]
![PNP问题][4]
![参考点][5]
![求解beta][6]


# 推导
> 见页面 [推导](https://www.zybuluo.com/snuffles/note/1455004)


  [1]: https://wx1.sinaimg.cn/mw690/c7716318ly1g24mz3ytzbj20cy05xq3j.jpg
  [2]: https://wx3.sinaimg.cn/mw690/c7716318ly1g24kwqwl4vj20i00ao3zr.jpg
  [3]: https://wx4.sinaimg.cn/mw690/c7716318ly1g24mo613j1j20e708itaq.jpg
  [4]: https://wx2.sinaimg.cn/mw690/c7716318ly1g265jck9glj20go0biwg9.jpg
  [5]: https://wx3.sinaimg.cn/mw690/c7716318ly1g265jckm8fj20fi0arabd.jpg
  [6]: https://wx2.sinaimg.cn/mw690/c7716318ly1g265jckikkj20f709s0to.jpg