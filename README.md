更正：
我们知道FBR——>DVR的本质其实就是一个表象变换，FBR的基组为sinc函数，DVR基组为格点$`|x_l\rangle`$，因此
```math
|\Psi_i\rangle=\sum_n c_{ni}|\phi_n\rangle=\sum_{nl}c_{ni}|x_l\rangle\langle x_l|\phi_n\rangle=\sum_{nl}c_{ni}B_{ln}/\sqrt{w_l}|x_{l}\rangle=\sum_l \tilde{c}_{li}|x_l\rangle
```
因此对DVR波函数作图时，不能只对$`\tilde{c}_{li}`$作图，而要除以$`\sqrt{w_l}`$，因为
```math
\tilde{c}_{li}=\sqrt{w_l}\langle x_l|\Psi_{i}\rangle=\sqrt{w_l}\Psi_{i}(x_l)
```
我作图时并未进行此步操作，特此更正！
