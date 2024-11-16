更正：
我们知道FBR——>DVR的本质其实就是一个表象变换，FBR的基组为sinc函数，DVR基组为格点$`|x_l\rangle`$，因此
```math
|\Psi_i\rangle=\sum_n c_{nl}|\phi_n\rangle=\sum_{nl}c_{nl}|x_l\rangle\langle|x_l|\phi_n\rangle=\sum_{nl}B_{ln}|x_{l}\rangle=\sum_l c_{ln}^{\prime}|x_l\rangle
```
因此对DVR波函数作图时，不能只对$`c_{ln}^{\prime}`$作图，而要除以$`\sqrt{w_l}`$，因为
```math
c_{ln}^{\prime}=\langle x_l|\Psi_{n}\rangle=\sqrt{w_l}\Psi_{n}(x_l)
```
这一步用了高斯积分，我作图时并未进行此步操作，特此更正！
