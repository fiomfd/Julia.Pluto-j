### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5ffc5ec0-f91f-11ec-0552-37ef5f25102d
begin
    using Symbolics
	using ForwardDiff
	using PlutoUI
    using LaTeXStrings
	using Formatting
	using Plots
	using Colors 
	using ImageFiltering
	using Images
	using ImageMagick
	using TestImages
	using ImageView
	using ImageTransformations
	using LinearAlgebra
	using LowRankApprox
	using Wavelets
	using DSP
	using Primes
    using DataFrames
	using ImplicitPlots
end

# ╔═╡ f3034480-e17d-4a7b-bdf3-cfb1d55d7bb9
md"""
## 2025年4月28日(月)
# #３ 計算論的思考とデータ
### 千原浩之

## 8:30-9:30 Zoom 講義
### I. 計算論的思考 (computational thinking)
### II. 線形代数
### III. 様々なデータとデータ処理

## 9:40-9:55 WebClass 試験

"""

# ╔═╡ 0a085483-72c0-4ec5-a677-9cf09f05cc54


# ╔═╡ 2aa2a4f4-65d4-43e4-b8fb-13bcca1b84b9
md"""
#### 計算論的思考 (computational thinking) とは何か？

- 数学および計算機を積極的に援用して問題の可視化や解決を図るための定式化を行うことである.


- 現在では世界中の大学で古典的数学とプログラミング言語の両方に長けた教員による "computational thinking" という題目の講義科目が提供されるようになっており、現代の必須の教養である.


- [MIT 18.C25 Julia: Solving Real-World Problems with Computation](https://computationalthinking.mit.edu/Fall24/) が世界で最も優れた computational thinking の講義の実践例として有名であり世界標準になっているが, 世界で最も優れた教育工学の実践例としても有名である.


- 日本国では高等学校情報Iの教科書の序文に登場する言葉でしかなく, あまり知られていない言葉である. 残念ながら古典的数学とプログラミング言語の両方に長けた教員があまりいないので, このような講義を提供することができない. また日本語社会では数学とは一切関係のない「計算論的思考の偽物」が出回っているので要注意である.
"""

# ╔═╡ 7d7e3384-cce9-4a97-bd78-604932418d56


# ╔═╡ 38a19ae1-07b9-414c-999c-c74d2b8df7d5
md"""
#### 微分可能性とは何かを計算論的思考する
"""

# ╔═╡ c2b971bd-5dc2-4f84-80d2-dd5307a96961
md"""
#### 1. 関数のグラフと直線の関係

区間$I$上の1変数関数 $f(x)$ のグラフ $\Gamma(f):=\{(x,f(x)) : x \in I\}$ 上の点 $(a,f(a))$, $a \in I$ を通る任意の傾き $\beta$ の直線 $y=g(x):=f(a)+\beta(x-a)$ を考えて, $f(x)$ と $g(x)$ のグラフの点 $(a,f(a))$ における近づき方 $f(x)-g(x)$ を観察する. 
"""

# ╔═╡ 763ca8d2-4975-4599-9aa2-247d2775c0cc


# ╔═╡ 84e0ec88-3824-4ad8-a38c-c01d9cacc7d7
md"""
#### ２. 関数のグラフに直線が接するとはどういうことか
高校数学では「接する」という言葉自体は高校数学IIの
- 円と直線の関係
- 微分可能性
に登場するが、それがどういうことなのかは何も説明がない。 実際には微分可能性の定義を論ずると「接する」ということの定義が見えてくる. 

ここでは原点を中心とする単位円と直線の関係を観察してみよう。 高校数学IIでは $x^2+y^2=1$ と方程式 $y=ax+b$ で表される直線が接するとは, 交点がただ一つの場合と定義されていてるだけであり, 交点付近での円弧と直線の振る舞い方についての説明はない. このような交点は $a^2+1=b^2$ が成立する場合にのみ起こり, 直線の方程式と交点のは次の2通りになる:

$\begin{aligned}
  y=ax+\sqrt{1+a^2}, 
& \quad
  \text{P}\left(-\frac{a}{\sqrt{1+a^2}},\frac{1}{\sqrt{1+a^2}}\right),
\\
  y=ax-\sqrt{1+a^2}, 
& \quad
  \text{Q}\left(\frac{a}{\sqrt{1+a^2}},-\frac{1}{\sqrt{1+a^2}}\right).
\end{aligned}$

P と Q を結ぶ直線の方程式は $y=-x/a$, 但し $a=0$ のときは $x=0$ で原点 $(0,0)$ を通り, 
$y=ax+b$ と直交することに注意する. 

 $a=\tan\theta$, $\theta\in(-\pi/2,\pi/2)$ とし, $\theta$ をパラメータとして図示する.

θ = $(@bind θ1 Slider(-89:1:89, show_value=true, default=30)) 度
"""

# ╔═╡ 5b8d6598-f602-4408-b948-84cec2bfc52a
begin
	f100(x,y) = x^2 + y^2 - 1
	P=zeros(2)
	Q=zeros(2)
	P[1]=-tan(θ1*pi/180)/sqrt(1+tan(θ1*pi/180)^2)
	Q[1]=1/sqrt(1+tan(θ1*pi/180)^2)
	P[2]=tan(θ1*pi/180)/sqrt(1+tan(θ1*pi/180)^2)
	Q[2]=-1/sqrt(1+tan(θ1*pi/180)^2)
	implicit_plot(f100;
				  grid=false,
				  xlim=(-2,2),
				  xlabel=L"x",
				  xticks = ([-2 -1 0 1 2;], [-2,-1,0,1,2]),
				  ylabel=L"y",
				  yticks = ([-2 -1 0 1 2;], [-2,-1,0,1,2]),
				  legend=false, 
				  linewidth=2)
	plot!([-2,2],[-2*tan(θ1*pi/180)+sqrt(1+tan(θ1*pi/180)^2),
				  2*tan(θ1*pi/180)+sqrt(1+tan(θ1*pi/180)^2)],
	 linecolor=:blue,
	 linewidth=2,
	 legend=:false)
	plot!([-2,2],[-2*tan(θ1*pi/180)-sqrt(1+tan(θ1*pi/180)^2),
				  2*tan(θ1*pi/180)-sqrt(1+tan(θ1*pi/180)^2)],
	 linecolor=:green,
	 linewidth=2,
	 legend=:false)
	plot!([-1.5*tan(θ1*pi/180)/sqrt(1+tan(θ1*pi/180)^2),
		  1.5*tan(θ1*pi/180)/sqrt(1+tan(θ1*pi/180)^2)],
		 [1.5/sqrt(1+tan(θ1*pi/180)^2),
		  -1.5/sqrt(1+tan(θ1*pi/180)^2)],
	 linecolor=:orange,
	 linewidth=1,
	 legend=:false)
	scatter!(P,Q, markersize = 5,label = false, markercolor=:magenta)
	scatter!((0,0), markersize = 5,label = false, markercolor=:red)
end

# ╔═╡ 05c00c05-c990-4632-ba43-4f6e14fad19f


# ╔═╡ 4c26a3e2-7ad7-4ca5-a5bb-1cf747f1a11a
md"""
#### 3. 予備的考察
p = $(@bind r Slider(0.1:0.1:10, show_value=true, default=2)) 
"""

# ╔═╡ c0233b8e-9db1-4a3a-8080-9f4bf7f3804b


# ╔═╡ 76cb2ec2-fc76-4ad5-8d57-777395755f7f
md"""
#### 4. $f(x)=x^2$  
$g(x):=f(a)+\beta(x-a)={\color{orange}{a^2+\beta(x-a)}}$

とすると

$
\begin{aligned}
  f(x)-g(x)
& =
  x^2-a^2-\beta(x-a)
\\
& =(x+a)(x-a)-\beta(x-a)
\\
& =
  (x+a-\beta)(x-a)
\\
& =
  (2a-\beta)(x-a)+(x-a)^2
\end{aligned}$

となる. 特に $\beta=2a$ のときの $g(x)$ を $h(x)={\color{green}{a^2+2a(x-a)}}$ とすると

$R(x):=
f(x)-h(x)=(x-a)^2$

となることに注意する. $\beta=2a$ のときのみ1次式 $x-a$ よりも速く $0$ に近づく.

a = $(@bind a1 Slider(-4:0.1:4, show_value=true, default=2)) 
β = $(@bind b1 Slider(-40:0.1:40, show_value=true, default=-1))
"""

# ╔═╡ 21c924bb-23de-42ac-805b-a19744a8e105


# ╔═╡ cfb114a5-e5c4-429d-bf6e-9acf7ea4b31d
md"""
#### 4. $f(x)=\lvert{x}\rvert$, $a=0$  
$g(x):=f(0)+\beta(x-0)={\color{orange}{\beta{x}}}$

とすると

$f(x)-g(x)
=
\lvert{x}\rvert-\beta{x}
=
\begin{cases}
-(\beta-1)x &\ (x \geqq0),
\\
-(\beta+1)x &\ (x <0).
\end{cases}$

任意の $\beta$ に対して $\beta\pm1$ の少なくとも一方は $0$ にならないので, $x-a$ よりも速く $0$ に近づくことはない. 

β = $(@bind b2 Slider(-40:0.1:40, show_value=true, default=0.3))
"""

# ╔═╡ 728a353b-e2fe-460d-8dcc-0f4cadbbd033


# ╔═╡ 2db362e9-e32a-4ba2-8810-3bdef318af41
md"""
#### 5. $f(x)=x^3-x$  
$g(x):=f(a)+\beta(x-a)={\color{orange}{a^3-a+\beta(x-a)}}$

とすると

$
\begin{aligned}
  f(x)-g(x)
& =
  x^3-x-a^3+a-\beta(x-a)
\\
& =(x^2+ax+a^2-1-\beta)(x-a)
\\
& =
  (3a^2-1-\beta)(x-a)+(x+2a)(x-a)^2
\end{aligned}$

となる. 特に $\beta=3a^2-1$ のときの $g(x)$ を $h(x)={\color{green}{a^3-a+(3a^2-1)(x-a)}}$ とすると

$R(x):=
f(x)-h(x)=(x+2a)(x-a)^2$

となることに注意する. $\beta=3a^2-1$ のときのみ1次式 $x-a$ よりも速く $0$ に近づく. 

a = $(@bind a3 Slider(-4:0.1:4, show_value=true, default=1.5)) 
β = $(@bind b3 Slider(-40:0.1:40, show_value=true, default=-1))
"""

# ╔═╡ b61b5be5-2b6e-4450-a3f1-6abe61bf9974


# ╔═╡ 2dc43448-3cf4-475f-8b61-d16c3594d9e3
md""" 
#### 6. 微分可能性の定義

1変数関数 $f(x)$ は $x=a$ の近くの区間 $(a-\delta,a+\delta)$ で定義されているものとする. 

関数$f(x)$ は $x=a$ で微分可能であるとは, 
- ある実数 $\alpha$ が存在して
- ある $x=a$ の近くで定義されている関数 $R(x)$ が存在して

$f(x)=f(a)+\alpha(x-a)+R(x), 
\qquad
\frac{R(x)}{x-a} \rightarrow 0 \quad (x \rightarrow a)$

が $x=a$ の近くで成立することと定義する. 

このとき, 実数 $\alpha$ を $f(x)$ の $x=a$ における微分係数とよび, 

$\frac{df}{dx}(a), \quad f^\prime(a)$

のように表す. また, 点 $\bigl(a,f(a)\bigr)$ を通って傾き $\alpha=f^\prime(a)$ の直線 

$y
=f(a)+f^\prime(a)(x-a)$

を $f$ のグラフ $\Gamma(f)$ の点 $\bigl(a,f(a)\bigr)$ における接線という. 

標語的に言えば、$f(x)$ が $x=a$ で微分可能であるとは, 
$f$ のグラフ $\Gamma(f)$ の点 $\bigl(a,f(a)\bigr)$ における接線がただ一つ存在することである. 
"""

# ╔═╡ 06191b8b-ddcc-4a2b-8228-0680aa8ad7c5
md""" 
#### 7. 定理 (微分可能性の特徴付け)
関数$f(x)$ は $x=a$ で微分可能であることは次と同値である. 

ある実数 $\alpha$ が存在して次が成立する:

$\frac{f(x)-f(a)}{x-a} 
\rightarrow \alpha
\quad
(x \rightarrow a).$ 
"""

# ╔═╡ b5e407b0-f3b2-4c9c-9274-9f4a69660110
md"""
#### 8. 定理の証明

必要性: 
関数$f(x)$ は $x=a$ で微分可能であるとせよ。すなわち、ある実数 $\alpha$ と関数 $R(x)$ が存在して

$f(x)=f(a)+\alpha(x-a)+R(x), 
\qquad
\dfrac{R(x)}{x-a} \rightarrow 0 \quad (x \rightarrow a)$

が $x=a$ の近くで成立するとせよ. このとき, 次が成立する: 

$\frac{f(x)-f(a)}{x-a} - \alpha = \frac{R(x)}{x-a} \rightarrow 0 
\quad (x \rightarrow a).$


十分性: 
ある実数 $\alpha$ が存在して, 次が成立するとせよ:  

$\frac{f(x)-f(a)}{x-a} \rightarrow \alpha  
\quad (x \rightarrow a).$

ここで

$R(x):=
f(x)-f(a)-\alpha(x-a)$

と定義すると, 次が成立する: 

$f(x)=f(a)+\alpha(x-a)+R(x),
\qquad
\dfrac{R(x)}{x-a}
=
\frac{f(x)-f(a)}{x-a}-\alpha \rightarrow 0
\quad
(x \rightarrow a).$
"""

# ╔═╡ f64936c9-a313-4c4b-a0b1-5e0df1b9e508


# ╔═╡ dbc102fd-f90e-4083-9155-173a0fd090a5
function newton1D(f, x0)
	
	f′(x) = ForwardDiff.derivative(f, x)   # \prime<TAB>
	
	x0 = 37.0  # starting point
	sequence = [x0]
	
	x = x0
	
	for i in 1:10
		x -= f(x) / f′(x)
	end
	
	return x
	
end

# ╔═╡ 21fda1bd-15c8-40e9-bb4f-ebbbd6c102ef
function newton2D_step(T, x)
	
	J = ForwardDiff.jacobian(T, x)   # should use StaticVectors
	
	δ = J \ T(x)   # J^(-1) * T(x)
	
	return x - δ
end

# ╔═╡ 7bb3421a-89ab-4566-aafc-c07350243d9c
"Looks for x such that f(x) = y, i.e. f(x) - y = 0"
function inverse(f, y, x0=[0, 0])
	return newton2D(x -> f(x) - y, x0)
end

# ╔═╡ df9b6dd3-9e44-4479-8049-233c774d345d
straight(x0, y0, x, m) = y0 + m * (x - x0)

# ╔═╡ 0ce69d78-5b11-4012-91e2-cbb9bde9e22f
function standard_Newton(f, n, x_range, x0, ymin=-10, ymax=10)
    
    f′ = x -> ForwardDiff.derivative(f, x)


	p = plot(f, x_range, lw=3, ylim=(ymin, ymax), legend=:false, size=(400, 300))

	hline!([0.0], c="magenta", lw=3, ls=:dash)
	scatter!([x0], [0], c="green", ann=(x0, -5, L"x_0", 10))

	for i in 1:n

		plot!([x0, x0], [0, f(x0)], c=:gray, alpha=0.5)
		scatter!([x0], [f(x0)], c=:red)
		m = f′(x0)

		plot!(x_range, [straight(x0, f(x0), x, m) for x in x_range], 
			  c=:blue, alpha=0.5, ls=:dash, lw=2)

		x1 = x0 - f(x0) / m

		scatter!([x1], [0], c="green", ann=(x1, -5, L"x_%$i", 10))
		
		x0 = x1

	end

	p |> as_svg


end

# ╔═╡ 4d7243c8-126d-40ac-a579-e53406a102d9
md"""
#### 微分積分に現れるその他の可視化
"""

# ╔═╡ 6efa7afe-b600-4ce4-81cb-43aa7348d8c2
md"""
#### 9. 1変数関数の零点を求めるニュートン法
関数 $f(x)$ は各点で微分可能とする. $f(x)$ の零点を求めることを考える. 
数列 $\{x_n\}$ を 

$x_{n+1}:=x_n-\frac{f(x_n)}{f^\prime(x_n)}, \quad n=0,1,2,\dotsc$

と定義すると, ある実数 $c$ が存在して

$x_n \rightarrow c \quad (n \rightarrow \infty), \qquad f(c)=0$

となることがある. このようにして $f(c)=0$ となる $c$ を求める方法をニュートン法という. 
"""

# ╔═╡ 4edf5113-e3a3-4f00-9fbc-1aebe02dfed1
md"""
#### $f(x)=x^2-9$
n = $(@bind n2 Slider(0:10, show_value=true, default=0))
x₀ = $(@bind x02 Slider(-4:10, show_value=true, default=6))
"""

# ╔═╡ baa5dd97-3152-469f-a45e-d427c786dcfd
let
	f(x) = x^2 - 9

	standard_Newton(f, n2, -4:0.01:10, x02, -15, 70)
end

# ╔═╡ 407adc66-f13c-4de1-8ad8-f59bf6b66063


# ╔═╡ 71e8ddf6-8c83-46ad-9222-3f1a9f57f2bf
md"""
#### 10. $\sin{x}$ のテイラー展開のコンパクト一様収束

実数直線上の関数 $\sin{x}$ は次のべき級数で表されることが知られている:

$\sin{x}=\sum_{k=0}^\infty\frac{(-1)^kx^{2k+1}}{(2k+1)!}=x-\frac{x^3}{6}+\frac{x^5}{120}-\dotsb, \quad x\in\mathbb{R}.$

右辺は $\sin{x}$ の $x=0$ を中心とするテイラー展開とよばれる. 右辺のべき級数は任意の有界閉区間上で一様収束して $\sin{x}$ に一致することが知られている. すなわち任意の $R>0$ に対して

$\max_{x\in[-R,R]}
\left\lvert
\sin{x}
-
\sum_{k=0}^K
\frac{(-1)^kx^{2k+1}}{(2k+1)!}
\right\rvert
\rightarrow 0 \quad (K\rightarrow\infty)$

が成立する. これを観察しよう. 


"""

# ╔═╡ 352fa6e4-1fc2-48ba-a62f-0e60ca0a6304
begin
    y = range(-3*pi, 3*pi, length = 301);
	f4=sin.(y);
	K=15;
	s3=zeros(K,301);
	for k=1:K
	    for l=1:301
	        s3[k,l]=(-1)^(k-1)*y[l]^(2*k-1)/factorial(big(2*k-1));
        end
	end
	S4=zeros(K,301);
	S4[1,:]=s3[1,:];
	for l=1:301
	    for k=2:K
		    S4[k,l]=S4[k-1,l]+s3[k,l];
	    end
	end
end

# ╔═╡ 76fa77f0-d82e-40ac-b826-50b8c781d030
md"""
K = $(@bind l4 Slider(1:K, show_value=true))
"""

# ╔═╡ b742bed1-9ae5-44f7-8b32-6609cdc3741b
begin
	h4=S4[l4,:];
    plot([f4,h4],grid=false,linewidth=2,ylim=(-1.2,1.2),
		title="sin(x) and its Taylor series",
		xticks = ([0 50 100 151 201 251 301;], ["-3π","-2π","-π","0","π","2π","3π"]),
		xlabel=L"x",
		label=[L"\sin{x}" "Taylor"],
		legend=:topleft,legendfont=font(8))
end

# ╔═╡ 6facfa26-5ae1-484f-a4f0-e8846ee8aec4


# ╔═╡ 50841743-af5c-4cea-9e61-6514ee36ec72
md"""
#### 11. 1変数連続関数のリーマン和の収束
有界閉区間 $[a,b]$ 上の実数値連続関数 $f(x)$ が与えられているとする. $[a,b]$ を$n$等分し、分点を

$a_k=a+\dfrac{k(b-a)}{n}, \quad k=0,1,\dotsc,n$

とすると, $a_0=a$, $a_n=b$ である. 各小区間 $[a_{k-1},a_k]$ の真ん中の点の $f$ の値をとって小区間の幅をかけて $k$ について和をとったリーマン和

$
\begin{aligned}
  R_n
& :=
  \sum_{k=1}^n
  f\left(\dfrac{a_{k-1}+a_k}{2}\right)
  \cdot
  (a_k-a_{k-1}) 
\\
& =
  \sum_{k=1}^n
  f\left(a+\dfrac{(2k-1)(b-a)}{2n}\right)
  \cdot
  \frac{b-a}{n}
\end{aligned}$

を考える. これは各小区間ごとの長方形の符号付き面積を足し合わせたものである. 

リーマン和 $R_n$ は $n \rightarrow \infty$ のとき, ある実数 $\alpha$ に収束することが知られている:

$R_n \rightarrow \alpha \quad (n \rightarrow \infty).$

実数 $\alpha$ を $f(x)$ の $[a,b]$ 上のリーマン積分とよび

$\int_a^bf(x)dx$

と表す。この例を観察しよう。
"""

# ╔═╡ bd09658a-e155-42ef-bebe-5088f55f177f
md"""
n= $(@bind n5 Slider(1:300, show_value=true))
"""

# ╔═╡ d9b4f11f-3687-4872-8237-bcfb396e9d88
begin
	x5 = range(0, 1, length = 501);
	f5=sin.(2*pi*x5)+x5+1.3*ones(501);
	
	g5=zeros(501);
	for j=1:501
		g5[j]=f5[Int64(floor(floor(n5*j/501)*501/2/n5)+ceil(ceil(n5*j/501)*501/2/n5))];
	end
end

# ╔═╡ 7fbec145-da40-4636-875d-c5124c6da60b
begin
    plot(f5,
		title="Convergence of Riemann sum",
	    grid=false,
		ylim=(0,2.7),
	    xticks=([1 501;],[L"a",L"b"]),
	    yticks=([0;],[0]),
		xlabel=(L"x"),
		ylabel=(L"y"), 
		edgecolor=:false, 
	    legend=false,
	    lw=3,
	    color=:magenta)
	bar!(g5,lw=0,color=:blue)
end

# ╔═╡ 14b5204e-f53d-4f5c-aaa1-070dde6eed35


# ╔═╡ b469f032-77fa-48a3-976d-79de94ddea4a
md"""
#### 12. ニュートン法と勾配降下法
現実のAIの多くはより下位の概念である深層学習のある種の発展程度で説明される概念である. AIの設計は既存のデータによる学習によるが、学習とは変数の個数が1億個くらいの関数の極小点の探索に過ぎない. ここでは極小点の代表的な数値的探索法であるニュートン法の2変数関数の例で観察する. 
"""

# ╔═╡ 4843d19a-ef5d-47d5-8c23-8dedeb069641
begin
	N6=100;
	x6 = range(-1.7, length=27, stop=0.9);
	y6 = range(-1.3, length=27, stop=1.3);
	f6=zeros(length(x6),length(y6));
	for k=1:27
		for l=1:27
			f6[l,k]=2*x6[k]^3+x6[k]*y6[l]^2+5*x6[k]^2+y6[l]^2;
		end
	end
end

# ╔═╡ 5b538799-ea23-46cf-9f40-34e902cfc91d
md"""
 learning rate $s$ = $(@bind s6 Slider(0.01:0.01:0.2, show_value=true, default=0.02))

 $x_0$ = $(@bind a6 Slider(-1:0.1:1, show_value=true, default=0.8)) 

 $y_0$ = $(@bind b6 Slider(-1:0.1:1, show_value=true, default=0.8))

 steps $n$ = $(@bind n6 Slider(0:1:N6, show_value=true, default=0))
"""

# ╔═╡ feaf5998-8014-4ff6-b600-70af413ebbba
begin
	Y4=zeros(N6+1,3);
	Y4[1,1]=a6;
	Y4[1,2]=b6;
	Y4[1,3]=2*a6^3+a6*b6^2+5*a6^2+b6^2;
	for n=1:N6
		D6=12*Y4[n,1]^2-2*Y4[n,2]^2+22*Y4[n,1]+10;
		Y4[n+1,1]=(12*Y4[n,1]^3-2*Y4[n,1]*Y4[n,2]^2+21*Y4[n,1]^2+Y4[n,2]^2+9*Y4[n,1])/D6;
		Y4[n+1,2]=Y4[n,2]*(12*Y4[n,1]^2-Y4[n,2]^2+17*Y4[n,1]+5)/D6;
		Y4[n+1,3]=2*Y4[n+1,1]^3+Y4[n+1,1]*Y4[n+1,2]^2+5*Y4[n+1,1]^2+Y4[n+1,2]^2;
	end

	Y5=zeros(N6+1,3);
	Y5[1,1]=a6;
	Y5[1,2]=b6;
	Y5[1,3]=2*a6^3+a6*b6^2+5*a6^2+b6^2;
	for n=1:N6
		Y5[n+1,1]=Y5[n,1]-s6*(6*Y5[n,1]^2+Y5[n,2]^2+10*Y5[n,1]);
		Y5[n+1,2]=Y5[n,2]-s6*(2*Y5[n,1]*Y5[n,2]+2*Y5[n,2]);
		Y5[n+1,3]=2*Y5[n+1,1]^3+Y5[n+1,1]*Y5[n+1,2]^2+5*Y5[n+1,1]^2+Y5[n+1,2]^2;
	end
end

# ╔═╡ 2196331a-3546-47f8-a423-5c49441e2bbb
begin
   surface(x6,y6,f6,
		    size=(600,600),
	        grid=false,
		    xlim=(-1.7,0.9),
		    ylim=(-1.3,1.3),
		    zlim=(0,9),
		    alpha=0.8, 
	        title="Newton's Method vs Gradient Descent",
		    xticks=([-1 0 1;],[-1,0,1]),
		    yticks=([-1 0 1;],[-1,0,1]),
		    zticks=([0 5;],[0,5]),
		    xlabel="\$x\$",
		    ylabel="\$y\$",
	        colormap=:jet,
	        colorbar=false,
		    camera = (-20, 30))
	scatter!([0],[0],[0],color=:orange,label="Minimal Point",markersize=7)
	scatter!([Y4[1:n6+1,1]],[Y4[1:n6+1,2]],[Y4[1:n6+1,3]],color=:yellow,label="Newton's Method",markersize=4)
	scatter!([Y5[1:n6+1,1]],[Y5[1:n6+1,2]],[Y5[1:n6+1,3]],color=:magenta,label="Gradient Descent",markersize=4)
	scatter!([Y4[1:1,1]],[Y4[1:1,2]],[Y4[1:1,3]],color=:red,label="Starting Point",markersize=6)
end

# ╔═╡ d2d1002e-12d0-4470-8681-731003de6c6d


# ╔═╡ 431d018d-c008-47f3-9d30-7e9da54d579a


# ╔═╡ 8d9148e7-4ce1-4da7-9eb3-43df098d5ad5
begin
	p0=zeros(2)
	q0=zeros(2)
	p0[1]=1
	p0[2]=101
	q0[1]=0
	q0[2]=1
    t = range(0, 1.2, length = 121)
	supersub=zeros(121)
	linear=zeros(121)
	for k=1:121
		supersub[k]=t[k]^r
		linear[k]=t[k]
	end
end

# ╔═╡ 769adba4-34e2-4e63-9099-2df3f6dbb529
begin
    plot([supersub,linear],grid=false,
		linewidth=2,ylim=(-0.1,1.5),
        title=L"s=t^p, s=t",
        xticks = ([0 101;], [0,1]),
		xlabel=L"t",
		yticks = ([0 1;], [0,1]),
		ylabel=L"s",
	    label = [L"s=t^p" L"s=t"],
	    legend = :topleft)
	plot!([0, 121],[0, 0], linestyle = :dash, linewidth=1,legend=:false)
	scatter!(p0,q0,markersize = 5,label = false)
end

# ╔═╡ ab1c3d69-0499-4fad-b2a8-5360ed7339e1
begin
	p1=zeros(1)
	q1=zeros(1)
    p1[1]=501+100*a1
	q1[1]=a1^2
	p3=zeros(1)
	q3=zeros(1)
    p3[1]=501+250*a3
	q3[1]=a3^3-a3
	x1 = range(-5, 5, length = 1001)
	x3 = range(-2, 2, length = 1001)
	f1=zeros(1001)
	g1=zeros(1001)
	h1=zeros(1001)
	f2=zeros(1001)
	g2=zeros(1001)
	f3=zeros(1001)
	g3=zeros(1001)
	h3=zeros(1001)
	for k=1:1001
		f1[k]=x1[k]^2
	    g1[k]=b1*(x1[k]-a1)+a1^2
		h1[k]=a1^2+2*a1*(x1[k]-a1)
		f2[k]=abs(x1[k])
		g2[k]=b2*x1[k]
		f3[k]=x3[k]^3-x3[k]
	    g3[k]=b3*(x3[k]-a3)+a3^3-a3
		h3[k]=a3^3-a3+(3*a3^2-1)*(x3[k]-a3)
	end
end

# ╔═╡ 8ed1312b-f5a3-4220-a865-9a12ed09a691
begin
    plot([f1,g1,h1],grid=false,
		linewidth=2,ylim=(-2,22),
        title=L"y=x^2,  y=a^2+\beta(x-a),  y=a^2+2a(x-a)",
        xticks = :false,
		xlabel=L"x",
	    label = [L"y=x^2" L"y=a^2+\beta(x-a)" L"y=a^2+2a(x-a)"],
		ylabel=L"y",
		yticks=:false,
	    legend = :top)
	scatter!(p1,q1,markersize = 5,label = false)
end

# ╔═╡ 3c696822-5b77-4512-83f0-c524371cde5d
begin
    plot([f2,g2],grid=false,
		linewidth=2,ylim=(-1,4),
        title=L"y=|x|, y=\beta{x}",
        xticks = :false,
		xlabel=L"x",
	    label = [L"y=|x|" L"y=\beta{x}"],
		ylabel=L"y",
		yticks=:false,
	    legend = :top)
	scatter!((501,0),markersize = 5,label = false)
end

# ╔═╡ 6d2d8054-56e3-4499-a8ea-cfcf291361e8
begin
    plot([f3,g3,h3],grid=false,
		linewidth=2,ylim=(-7,7),
        title=L"y=x^3-x,  y=a^3-a+\beta(x-a),  y=a^3-a+(3a^2-1)(x-a)",
        xticks = :false,
		xlabel=L"x",
	    label = [L"y=x^3-x" L"y=a^3-a+\beta(x-a)" L"y=a^3-a+(3a^2-1)(x-a)"],
		ylabel=L"y",
		yticks=:false,
	    legend = :top)
	scatter!(p3,q3,markersize = 5,label = false)
end

# ╔═╡ 0e933e40-e613-4eda-bc8c-666a7d208c02
md"""
### II. 線形代数

共通教育科目 線形代数学I-II 01組 は数ベクトル空間に限定した線形代数を解説する. 
- 平面上の線形代数 (30年前までの高校数学IIBの話題)
- 連立1方程式の解の存在と解の全体の構造
- 行列式と計算法
- 数ベクトル空間
- 数ベクトル空間の間の線形写像
- 線形代数学の基本定理
- 正方行列の対角化
- 一般の型の行列の特異値分解とその応用 (線形写像の構造・画像の近似・線形回帰など)
- 線形写像の例: 有限離散フーリエ解析、離散ウェーブレットと多重解像度分解

[文部科学省：高等学校 行列入門](https://www.mext.go.jp/a_menu/shotou/new-cs/senseiouen/1394142_00001.html)
"""

# ╔═╡ b6af8639-aaf1-4148-9ddc-410d2a7b0495
md"""
##### 1. 行列と数ベクトルの定義

さて$m$と$n$を正整数とし$m{\times}n$個の(実)数 

$a_{ij}, \quad i=1,\dotsc,m, \quad j=1,\dotsc,n$ 

を行を $m$ 個で列が $n$ 個になるように並べて $[\quad]$ または $(\quad)$ で括ったもの

$A=[a_{ij}]
:=
\begin{bmatrix}
a_{11} & a_{12} & \dotsb & a_{1n}
\\
a_{21} & a_{22} & \dotsb & a_{2n}
\\
\vdots & \vdots &  & \vdots
\\
a_{m1} & a_{m2} & \dotsb & a_{mn}
\end{bmatrix}$

を $m{\times}n$ 行列という. 例えば

$\begin{bmatrix}
1 & 2 & 3
\\
4 & 5 & 6
\end{bmatrix},
\quad
\begin{bmatrix}
1 & 2 & 3
\\
4 & 5 & 6
\\
7 & 8 & 9
\end{bmatrix},
\quad
\begin{bmatrix}
1
\\
4
\\
7
\end{bmatrix},
\quad
\begin{bmatrix}
1 & 2 & 3
\end{bmatrix}$

は行列である. 特に列または行の数が1である場合

$\vec{a}
=
\begin{bmatrix}
a_1 \\ \vdots \\ a_m
\end{bmatrix},
\quad
\vec{b}
=
\begin{bmatrix}
b_1 & \dotsb & b_n
\end{bmatrix}$

はそれぞれ $m$ 次元列ベクトル, $n$ 次元行ベクトルとよぶ. 
"""

# ╔═╡ c15d7baf-54d5-49da-b144-3ce7255cc2f1
md"""
##### 2. 行列と数ベクトルの演算

- 行列 $A$ の転置行列 $A^T$: $m{\times}n$ 行列 $A=[a_{ij}]$ の転置行列 $A^T$ という $n{\times}m$ 行列を次のように定義する:
$A^T:=[a_{ji}], 
\quad 
\begin{bmatrix}
1 & 2 & 3
\\
4 & 5 & 6
\end{bmatrix}^T
=
\begin{bmatrix}
1 & 4
\\
2 & 5
\\
3 & 6
\end{bmatrix}.$

- 行列の和 $A+B$: $A=[a_{ij}]$ と $B=[b_{ij}]$ が同じ $m{\times}n$ 行列であるとき $A+B$ を次のように定義する:
$A+B:=[a_{ij}+b_{ij}],$
$\begin{bmatrix}
1 & 2 & 3
\\
4 & 5 & 6
\end{bmatrix}
+
\begin{bmatrix}
7 & 8 & 9
\\
10 & 11 & 12
\end{bmatrix}
=
\begin{bmatrix}
1+7 & 2+8 & 3+9
\\
4+10 & 5+11 & 6+12
\end{bmatrix}
=
\begin{bmatrix}
8 & 10 & 12
\\
14 & 16 & 18
\end{bmatrix}.$

- 行列のスカラー倍 $cA$: 数 $c$ と $A=[a_{ij}]$ に対して $A$ のスカラー倍 $cA$ を次のように定義する:
$cA:=[ca_{ij}],
\quad
5
\begin{bmatrix}
1 & 2 & 3
\\
4 & 5 & 6
\end{bmatrix}
=
\begin{bmatrix}
5 & 10 & 15
\\
20 & 25 & 30
\end{bmatrix}.$

- 行列と列ベクトルの積 $A\vec{x}$: $\vec{a}_1,\dotsc,\vec{a}_n$ を $m$ 次元列ベクトルとし, $m{\times}n$ 行列 $A$ を $A:=\begin{bmatrix}\vec{a}_1 & \dotsb & \vec{a}_n\end{bmatrix}$ とする. $n$ 次元列ベクトル $\vec{x}=\begin{bmatrix}x_1 \\ \vdots \\ x_n\end{bmatrix}$ に対して $m$ 次元列ベクトル $A\vec{x}$ を次のように定義する:
$A\vec{x}:=x_1\vec{a}_1+\dotsb+x_n\vec{a}_n,$
$\begin{bmatrix}
1 & 2
\\
3 & 4
\\
5 & 6
\end{bmatrix}
\begin{bmatrix}
100 \\ 10
\end{bmatrix}
=
100
\begin{bmatrix}
1
\\
3
\\
5
\end{bmatrix}
+
10
\begin{bmatrix}
2
\\
4
\\
6
\end{bmatrix}
=
\begin{bmatrix}
120
\\
340
\\
560
\end{bmatrix}.$

- 行列の積 $AB$: $m{\times}n$ 行列 $A$ と $n{\times}l$ 行列 $B:=\begin{bmatrix}\vec{b}_1 & \dotsb & \vec{b}_l\end{bmatrix}$ に対して $m{\times}l$ 行列 $AB$ を次のように定義する:
$AB:=\begin{bmatrix}A\vec{b}_1 & \dotsb & A\vec{b}_l\end{bmatrix},$
$\begin{bmatrix}
1 & 2
\\
3 & 4
\\
5 & 6
\end{bmatrix}
\begin{bmatrix}
10 & 10
\\ 
1 & 100
\end{bmatrix}
=
\begin{bmatrix}
12 & 210
\\
34 & 430
\\
56 & 650
\end{bmatrix}.$

"""


# ╔═╡ 3f17e529-3ca4-4302-a54b-30fe096893ac
md"""
##### 3. 平面ベクトルの加法

$\vec{a} = \begin{bmatrix} a_1 \\ a_2 \end{bmatrix}, \quad 
\vec{b} = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}, \quad
\vec{a} + \vec{b} = \begin{bmatrix} a_1 + b_1 \\ a_2 + b_2 \end{bmatrix}$

a1 = $(@bind aa1 Slider(-1.5:0.1:1.5, show_value=true, default=1.3))  
b1 = $(@bind bb1 Slider(-1.5:0.1:1.5, show_value=true, default=-1))  

a2 = $(@bind aa2 Slider(-1.5:0.1:1.5, show_value=true, default=0.3))  
b2 = $(@bind bb2 Slider(-1.5:0.1:1.5, show_value=true, default=1.2))  
"""

# ╔═╡ 157db428-8868-403b-9e3b-705f3df7612a
begin
    GR.setarrowsize(1.5)

    A1 = [0, aa1]; A2 = [0, aa2]
    B1 = [0, bb1]; B2 = [0, bb2]
    AB1 = [0, aa1+bb1]; AB2 = [0, aa2+bb2]
    C1 = [aa1, aa1+bb1, NaN, bb1, aa1+bb1]
    C2 = [aa2, aa2+bb2, NaN, bb2, aa2+bb2]

    plot(A1, A2,
        title = L"\vec{a} + \vec{b}",
        xlabel = L"x", ylabel = L"y",
        titlefontsize = 20,
        xlim = (-2.1, 2.1), ylim = (-2.1, 2.1),
        arrow = (:closed, 2.0),
        color = :blue, label = L"\vec{a}",
        grid = false, aspect_ratio = 1,
        legend = :bottomleft, legendfontsize = 14)

    plot!(B1, B2,
        arrow = (:closed, 2.0),
        color = :orange, label = L"\vec{b}")

    plot!(AB1, AB2,
        arrow = (:closed, 2.0),
        color = :green, label = L"\vec{a} + \vec{b}")

    plot!(C1, C2,
        line = (1, :dash, :purple),
        label = false)

    #annotate!(aa1, aa2, text(L"\vec{a}", :blue, :center, 10))
	#annotate!(bb1, bb2, text(L"\vec{b}", :orange, :center, 10))
    #annotate!(aa1+bb1, aa2+bb2, text(L"\vec{a}+\vec{b}", :green, :center, 10))
end


# ╔═╡ f95d5371-c27e-459f-9ca8-b36c4aa67177
md"""
##### 4. 平面上の線形写像
さて $A$ を $m{\times}n$ 行列, $\vec{x}$ と $\vec{y}$ を $n$ 次元列ベクトル, $c$ と $d$ を数とする. 
次の対応 

$\vec{x} \mapsto A\vec{x}$

は $n$ 次元列ベクトルの全体 $\mathbb{R}^n$ から $m$ 次元列ベクトルの全体 $\mathbb{R}^m$ への写像であり, 列ベクトルの和とスカラー倍を保存する: 

$A(c\vec{x}+d\vec{y}) = c(A\vec{x}) + d(A\vec{y}).$

このような性質をもつ写像は線形写像とよばれる. 以下では $2\times2$ 行列が定義する平面 $\mathbb{R}^2$ から平面 $\mathbb{R}^2$ への線形写像の例を2つ見てみよう.

$\vec{x} \mapsto \vec{s}:=A\vec{x},$

すなわち

$\begin{bmatrix} x \\ y \end{bmatrix}
\mapsto 
\begin{bmatrix} s \\ t \end{bmatrix}
:=
\begin{bmatrix} a & b \\ c & d \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix}
=
x
\begin{bmatrix} a \\ c \end{bmatrix}
+
y
\begin{bmatrix} b \\ d \end{bmatrix}$

とする.
"""

# ╔═╡ fbad5a13-c986-4d2a-936d-51757d7eac93
md"""
##### 5. 伸長
２つの定数 $\lambda\ne0$ かつ $\mu\ne0$ を導入する. 

$\begin{bmatrix} s \\ t \end{bmatrix}
:=
\begin{bmatrix} \lambda & 0 \\ 0 & \mu \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix}
=
\begin{bmatrix} \lambda{x} \\ \mu{y} \end{bmatrix}$

は $x$ 軸の尺度を $\lambda$ 倍に引き伸ばし, $y$ 軸の尺度を $\mu$ 倍に引き伸ばす写像である. 
$\lambda<0$ ならば $x$ 軸と $s$ 軸の向きは逆になる. 

原点を中心とする単位円の方程式 $x^2+y^2=1$ は楕円の方程式になる:

$\frac{s^2}{\lambda^2}+\frac{t^2}{\mu^2}=1$

λ = $(@bind lambda Slider(0.5:0.1:2, show_value=true, default=1))  
μ = $(@bind mu Slider(0.5:0.1:2, show_value=true, default=1))  

"""

# ╔═╡ 1299f123-8d64-423f-bc64-6410a7392cc2
begin
	f(s,t) = s^2/lambda^2 + t^2/mu^2 -4
	implicit_plot(f; xticks=false, yticks=false, legend=false, linewidth=3,xlabel = L"s", ylabel = L"t" )
end

# ╔═╡ 1c624c90-99ed-414d-a6e4-b2f6b16101c8
md"""
##### 6. 回転
原点を中心に反時計回りに角度 $\theta \in [0,2\pi)$ の回転を与える行列 $P(\theta)$ は

$P(\theta)
=
\begin{bmatrix}
\cos\theta & -\sin\theta
\\
\sin\theta & \cos\theta
\end{bmatrix}$

である. 第1列と第2列の関係は以下のようになり, 直交していることがわかる:

$\begin{bmatrix}
-\sin\theta
\\
\cos\theta
\end{bmatrix}
=
\begin{bmatrix}
\cos(\theta+\pi/2)
\\
\sin(\theta+\pi/2)
\end{bmatrix}.$

$st$-平面に $x$-軸と $y$-軸を書き込んでみる:

$\begin{bmatrix} s \\ t \end{bmatrix}
=
\begin{bmatrix}
\cos\theta & -\sin\theta
\\
\sin\theta & \cos\theta
\end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix}
=
x
\begin{bmatrix} \cos\theta \\ \sin\theta \end{bmatrix}
+
y
\begin{bmatrix} -\sin\theta \\ \cos\theta \end{bmatrix}.$

θ = $(@bind θ Slider(0:1:360, show_value=true, default=30))  

"""

# ╔═╡ 0969f38b-f1f6-4b81-b358-c87ab486beff
begin
    GR.setarrowsize(1.5)

    s1 = [-2, 2]; s2 = [0, 0]
    t1 = [0, 0]; t2 = [-2, 2]
    xx1 = [-2*cos(θ*pi/180), 2*cos(θ*pi/180)]; xx2 = [-2*sin(θ*pi/180), 2*sin(θ*pi/180)]
    yy1 = [2*sin(θ*pi/180), -2*sin(θ*pi/180)]; yy2 = [-2*cos(θ*pi/180), 2*cos(θ*pi/180)]
	
    plot(s1,s2,
        xlabel = L"s", ylabel = L"t",
		xlabelfontsize = 14,
		ylabelfontsize = 14,
        titlefontsize = 20,
        xlim = (-2.1, 2.1), ylim = (-2.1, 2.1),
		xticks=false, yticks=false,
        arrow = (:closed, 2.0),
        color = :blue, label = false,
        grid = false, aspect_ratio = 1,
        legend = :bottomleft, legendfontsize = 14)

    plot!(t1,t2,
        arrow = (:closed, 2.0),
        color = :blue, label = false)

    plot!(xx1,xx2,
        arrow = (:closed, 2.0),
        color = :green, label = L"x")

    plot!(yy1,yy2,
        arrow = (:closed, 2.0),
        color = :magenta, label = L"y")

end


# ╔═╡ 3e57202c-4c53-4282-a477-410a990bb318
md"""
##### 7. 鏡映
 $\theta \in (-\pi,\pi)$ とする. 次の線形写像

$\vec{v}=M(\theta)\vec{u},
\quad
\vec{u}
=
\begin{bmatrix} u_1 \\ u_2 \end{bmatrix},
\quad
\vec{v}
=
\begin{bmatrix} v_1 \\ v_2 \end{bmatrix},
\quad
M(\theta)
=
\begin{bmatrix}
\cos(2\theta) & \sin(2\theta)
\\
\sin(2\theta) & -\cos(2\theta)
\end{bmatrix}$

は直線 $y=x\tan\theta$ についての対称な点への写像である. 

u1 = $(@bind u1 Slider(-1.6:0.1:1.6, show_value=true, default=1.4))  
u2 = $(@bind u2 Slider(-1.6:0.1:1.6, show_value=true, default=1.4))

θ = $(@bind φ Slider(-90:1:90, show_value=true, default=20))  
"""

# ╔═╡ 659ab5e5-5b78-4496-a3d5-5f4b22396178
begin
	GR.setarrowsize(1.5)

    ax = [-3*sqrt(2)*cos(φ*pi/180), 3*sqrt(2)*cos(φ*pi/180)]
    ay = [-3*sqrt(2)*sin(φ*pi/180), 3*sqrt(2)*sin(φ*pi/180)]
	v1=u1*cos(φ*pi/90)+u2*sin(φ*pi/90)
	v2=u1*sin(φ*pi/90)-u2*cos(φ*pi/90)
    bx=[u1,v1]
	by=[u2,v2]
	
    plot(ax,ay,
        xlabel = L"x", ylabel = L"y",
		xlabelfontsize = 14,
		ylabelfontsize = 14,
        titlefontsize = 20,
        xlim = (-2.5, 2.5), ylim = (-2.5, 2.5),
		xticks=false, yticks=false,
        line = (:closed, 1.5),
        color = :green, label = false,
        grid = false, aspect_ratio = 1,
        legend = :bottomleft, legendfontsize = 12)
	plot!(bx,by,
        line = (:closed, 1.0),
        color = :orange, label = false)
	scatter!((u1,u2),markersize = 5,label=L"\vec{u}", color=:magenta)
	scatter!((v1,v2),markersize = 5,label=L"\vec{v}", color=:blue)
	scatter!((0,0),markersize = 5,label=L"\vec{0}", color=:red)
end

# ╔═╡ ca6d79ef-d5a4-4a55-b7f9-83d8682ddd1b
md"""
##### 8. 特異値分解

任意の $2\times2$ 行列 $A=\begin{bmatrix} a & b \\ c & d\end{bmatrix}$ に対して, ある $\theta, \varphi \in [0,2\pi]$ と $\lambda \geqq \mu \geqq 0$ が存在して

$A
=
\begin{bmatrix}
\cos\theta & -\sin\theta
\\
\sin\theta & \cos\theta
\end{bmatrix}
J
\begin{bmatrix}
\lambda & 0
\\
0 & \mu
\end{bmatrix}
\begin{bmatrix}
\cos\varphi & -\sin\varphi
\\
\sin\varphi & \cos\varphi
\end{bmatrix}
K$

が成立する. これを $A$ の特異値分解という. 特異値分解は一般の $m{\times}n$ 行列に一般化されるされる。 ここに

$J, K = 
\begin{bmatrix}
1 & 0 
\\
0 & 1
\end{bmatrix},
\begin{bmatrix}
0 & 1 
\\
1 & 0
\end{bmatrix}.$


行列 

$A
=
\frac{1}{8}
\begin{bmatrix} 
\sqrt{6}+4\sqrt{2} & \sqrt{6}-4\sqrt{2} 
\\ 
-4\sqrt{6}+\sqrt{2} & 4\sqrt{6}+\sqrt{2} 
\end{bmatrix}$ 

は次のように分解される:

$A
=
U \Sigma V^T
=
\begin{bmatrix} \cos(\pi/3) & \sin(\pi/3) \\ -\sin(\pi/3) & \cos(\pi/3) \end{bmatrix}
\begin{bmatrix} 2 & 0 \\ 0 & 1/2 \end{bmatrix}
\begin{bmatrix} \cos(\pi/4) & \sin(\pi/4) \\ -\sin(\pi/4) & \cos(\pi/4) \end{bmatrix}^T.$

さて $U$ は時計回りの $\pi/3$ の回転を与える行列であり, $V^T$ は時計回りの $\pi/4$ の回転を与える行列であり, $A$ は回転 ($V^T$), 伸長 ($\Sigma$), 回転 ($U$) の合成写像であることがわかる. 

$\begin{bmatrix} s \\ t \end{bmatrix}
=
A
\begin{bmatrix} x \\ y \end{bmatrix}
=
\frac{x}{8}
\begin{bmatrix} \sqrt{6}+4\sqrt{2} \\ -4\sqrt{6}+\sqrt{2} \end{bmatrix}
+
\frac{y}{8}
\begin{bmatrix} \sqrt{6}-4\sqrt{2} \\ 4\sqrt{6}+\sqrt{2}  \end{bmatrix}$

とすると, 単位円の方程式 $x^2+y^2=1$ は傾いた楕円の方程式 $49s^2+15\sqrt{3}st+19t^2=16$ に移る.  
"""

# ╔═╡ b06741a8-5ca6-4e15-b085-e3bbcddca93d
begin 
	svdexample=load("./image/svd.png");
end

# ╔═╡ 87dd79c0-c0e9-4021-879a-c18cfa6e9252
md"""
$\sigma_1=2,
\quad
\sigma_2=\frac{1}{2},$

$\vec{v}_1
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 \\ -1 
\end{bmatrix},
\quad
\vec{v}_2
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 \\ 1 
\end{bmatrix},$

$\vec{u}_1
=
\frac{1}{2}
\begin{bmatrix}
1 \\ -\sqrt{3} 
\end{bmatrix},
\quad
\vec{u}_2
=
\frac{1}{2}
\begin{bmatrix}
\sqrt{3} \\ 1 
\end{bmatrix}$

とおくと, 


$A
=
\sum_{j=1}^2\sigma_j\vec{u}_j\vec{v}_j^T
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1 \\ -\sqrt{3} 
\end{bmatrix}
\begin{bmatrix}
1 & -1 
\end{bmatrix}
+
\frac{1}{4\sqrt{2}}
\begin{bmatrix}
\sqrt{3} \\ 1 
\end{bmatrix}
\begin{bmatrix}
1 & 1 
\end{bmatrix},$

すなわち

$A
=
\frac{1}{8}
\begin{bmatrix} 
\sqrt{6}+4\sqrt{2} & \sqrt{6}-4\sqrt{2} 
\\ 
-4\sqrt{6}+\sqrt{2} & 4\sqrt{6}+\sqrt{2} 
\end{bmatrix}
=
\frac{1}{8}
\begin{bmatrix} 4\sqrt{2} & -4\sqrt{2} \\ -4\sqrt{6} & 4\sqrt{6} \end{bmatrix}
+
\frac{1}{8}
\begin{bmatrix} \sqrt{6} & \sqrt{6} \\ \sqrt{2} & \sqrt{2} \end{bmatrix}$

と表すことができる. 
"""

# ╔═╡ f68a62e7-f442-49d1-bd1d-b7268d55f037


# ╔═╡ 3e301f5d-cdb8-43c2-8562-a83416c270cb
md"""
### III. 様々なデータとデータ処理
"""

# ╔═╡ 173a1432-4ca9-4ce1-9de2-2ef441b1810d
begin
	xx=rand(1550:1950, 100)/10;
    XX=[ones(100) xx];
    yy=(xx.^2).*rand(160:400, 100)/100000;
	yy=round.(yy, digits=1);
    bb=XX\yy;
    tt=155:1:195;
    zz=[ones(length(tt)) tt]*bb;
	ZZ=[xx yy];
	for i=1:1
	end
end

# ╔═╡ 3cd15536-0761-44db-a8d1-96096d24e7c5
md"""
#### 1. データ行列
下記は100人の身長と体重のデータ一覧であるが, $100\times2$ の行列とみなすことができる. 一般に個体数 $N$ で各個体がそれぞれ $p$ 種類のデータをもっているとき, そのデータ一覧は $N{\times}p$ 行列として表すことができる. 
"""

# ╔═╡ 213f0010-d7b9-4584-8484-868c5cd9316d
begin
	df = DataFrame(height = xx, weight = yy)
end

# ╔═╡ f7a1f0f0-2cd1-412b-98e2-9b3bb71b5c43
md"""
身長を横軸に体重を縦軸にとった平面上に各個人のデータを点で表すことができる. これを散布図という. 線形代数の一般逆行列あるいは特異値分解を利用すると, 線形回帰(最小二乗法)による身長と体重の関係を表す直線を求めることができる.
"""

# ╔═╡ 295e359a-ae03-4d15-b6f9-ee7c98067eb9
begin
	scatter(xx,yy, 
		grid=false,
		label="sample", 
		xlabel="height [cm]", 
		ylabel="weight [kg]", 
        title="Linear Regression between Height & Weight")
    plot!(tt,zz,linewidth=2, label="regression",color="magenta")
end

# ╔═╡ 3315d5ae-444e-4958-a1d7-086eeaaad522


# ╔═╡ 7b45973d-bfce-4b9b-8d2c-e35b5ef421a5
md"""
#### 2. 白黒画像 (grayscale image)
さて, 白黒画像 (grayscale image) は各小正方形(ピクセル)に $0$ から $255$ までの整数, あるいは, 
それらを $255$ で割った $0$ から $1$ までの実数を成分とする行列である. 
$0$ は黒を表し $255$ は白を表す. 数値は色の強さを表し大きいほど強い. 例えば次の行列

$\begin{bmatrix}
0 & 15  & 30  & 45 & 60 & 75
\\ 
90 & 105 & 120  & 135  & 150 & 165
\\ 
180 & 195  & 210  & 225 & 240 & 255
\end{bmatrix}$

の表す白黒画像は以下のとおりである。
"""

# ╔═╡ 2a819b34-8fdc-47b9-9397-72e8f791a902
begin
    Gsample=[0 15 30 45 60 75; 
		     90 105 120 135 150 165;
		     180 195 210 225 240 255]/255;
	plot(Gray.(Gsample),
		xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false)
end

# ╔═╡ 6450807b-e3e9-43e6-bd64-50c385c077ce


# ╔═╡ fe18cdd8-c5d6-443f-b99d-66514a9393b8
md"""
#### 3. RGB画像

RGB画像は同じサイズの白黒画像のデータの3つの行列をデータとし, 各行列を三原色の赤と緑と青でそれぞれ着色して合わせたものである. プログラミング言語の乱数機能で生成した3つの$16\times16$行列とその3つ組が表すRGB画像を示す.

"""

# ╔═╡ 1de04ede-fa17-4497-a85d-2205df043d04
begin
# read image, resize, decompose}
	L1=4;
	RII=rand(2^L1, 2^L1);
	GII=rand(2^L1, 2^L1);
	BII=rand(2^L1, 2^L1);
	RRR=zeros(3,2^L1, 2^L1);
	GGG=zeros(3,2^L1, 2^L1);
	BBB=zeros(3,2^L1, 2^L1);
	RGBII=zeros(3,2^L1, 2^L1);
	RRR[1,:,:]=RII[:,:];
	GGG[2,:,:]=GII[:,:];
	BBB[3,:,:]=BII[:,:];
	RGBII[1,:,:]=RII[:,:];
	RGBII[2,:,:]=GII[:,:];
	RGBII[3,:,:]=BII[:,:];

QQQ1=plot(colorview(RGB,RRR),
        title="Grayscale Image (Red)",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
QQQ2=plot(colorview(RGB,GGG),
        title="Grayscale Image (Green)",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
QQQ3=plot(colorview(RGB,BBB),
        title="Grayscale Image (Blue)",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
QQQ4=plot(colorview(RGB,RGBII),
        title="RGB Image",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
plot(QQQ1,QQQ2,QQQ3,QQQ4,size=(1000,1000),layout=(2,2))
end

# ╔═╡ 0f18af6e-c2ba-43df-8299-2c42821c8c17


# ╔═╡ de50b16d-e218-4221-ba1a-a4539dc33973
md"""
#### 4. 白黒動画 (grayscale movie)
複数の白黒画像を連続的に流すと白黒動画になる。ここでは $5\times10$ の行列を乱数によって100個生成したものをアニメーションにしてみる.
"""

# ╔═╡ 02d137fc-4aac-4da4-8850-04a6529b7042
begin
# an example of three-dimensional data
Agray=rand(5,10,100);
anim1 = @animate for k=1:100
    plot(Gray.(Agray[:,:,k]),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false)
end
gif(anim1, "grayscale_animation.gif", fps = 3)	
end

# ╔═╡ 9f9c9544-7454-4277-84a4-7a2b46fdea2f


# ╔═╡ 8c938021-1a02-44f6-93b0-a474f89a60d4
md"""
#### 5. RGB動画 (RGB movie)
"""

# ╔═╡ 8c14a644-b7f3-45eb-bc0c-bc9fc28d65bd
begin
# an example of four-dimensional data
RR=rand(80:255, 5, 10, 100);
GG=rand(80:255, 5, 10, 100);
BB=rand(80:255, 5, 10, 100);
V=zeros(3,5,10,100);
VR=zeros(3,5,10,100); 
VG=zeros(3,5,10,100);
VB=zeros(3,5,10,100);
ZE=zeros(5,10,100);
V[1,:,:,:]=RR/255;
V[2,:,:,:]=GG/255;
V[3,:,:,:]=BB/255;
VR[1,:,:,:]=RR/255;
VR[2,:,:,:]=ZE;
VR[3,:,:,:]=ZE;
VG[1,:,:,:]=ZE;
VG[2,:,:,:]=GG/255;
VG[3,:,:,:]=ZE;
VB[1,:,:,:]=ZE;
VB[2,:,:,:]=ZE;
VB[3,:,:,:]=BB/255;
anim2 = @animate for k=1:100
    P1=plot(colorview(RGB,VR[:,:,:,k]),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false,
         title="Grayscale (Red)");
    P2=plot(colorview(RGB,VG[:,:,:,k]),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false,
         title="Grayscale (Green)");
    P3=plot(colorview(RGB,VB[:,:,:,k]),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false,
         title="Grayscale (Blue)");
    P4=plot(colorview(RGB,V[:,:,:,k]),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false,
         title="RGB");
   plot(P1,P2,P3,P4,
         layout=(2,2),
         xaxis=false, 
         xticks=false, 
         yaxis=false, 
         yticks=false, 
         grid=false)
end
gif(anim2, "rgb_animation.gif", fps = 3)	
end

# ╔═╡ 4d3d1122-2424-478f-894d-f7180e59fca2
begin
# read image, resize, decompose
	#I=load("./material/char_kway_teow.jpg");
	#X=imresize(I, ratio=1/4);
	I=load("./material/CityU.jpg");
	X=imresize(I, ratio=1/8);	
    (p,q)=size(X);
	L=7;
    A=channelview(X);
    R=Array{Float64}(A[1,:,:]);
	G=Array{Float64}(A[2,:,:]);
	B=Array{Float64}(A[3,:,:]);


	RU, RS, RV=psvd(R);
	GU, GS, GV=psvd(G);
	BU, BS, BV=psvd(B);
	
	rank=100;
	DR=zeros(p,q,rank);
	DG=zeros(p,q,rank);
	DB=zeros(p,q,rank);
    for r=1:rank
	    DR[:,:,r]=sum(RS[n]*RU[1:p,n]*(RV[1:q,n])' for n=1:r);
        DG[:,:,r]=sum(GS[n]*GU[1:p,n]*(GV[1:q,n])' for n=1:r);
        DB[:,:,r]=sum(BS[n]*BU[1:p,n]*(BV[1:q,n])' for n=1:r);
	end
end

# ╔═╡ 8ceef78f-5058-4e6a-b5ce-ea5b22f6a7fb
begin
    Y=zeros(3,p,q,rank+16);
	Z=zeros(p,q)
	for r=1:4
	    Y[1,:,:,r]=R;
        Y[2,:,:,r]=G;
        Y[3,:,:,r]=B;
	end
    for r=5:8
	    Y[1,:,:,r]=Z;
        Y[2,:,:,r]=Z;
        Y[3,:,:,r]=Z;
    end
	for r=9:rank+8
        Y[1,:,:,r]=DR[:,:,r-8];
        Y[2,:,:,r]=DG[:,:,r-8];
        Y[3,:,:,r]=DB[:,:,r-8];
	end
    for r=rank+9:rank+12
        Y[1,:,:,r]=Z;
        Y[2,:,:,r]=Z;
        Y[3,:,:,r]=Z;
	end
    for r=rank+13:rank+16
        Y[1,:,:,r]=R;
        Y[2,:,:,r]=G;
        Y[3,:,:,r]=B;
	end
end

# ╔═╡ 83b7d146-00a0-47b7-8988-288721afad10
begin
    W=zeros(3,p,q,rank);
    for r=1:rank
	    W[1,:,:,r]=DR[:,:,r];
        W[2,:,:,r]=DG[:,:,r];
        W[3,:,:,r]=DB[:,:,r];
	end
end

# ╔═╡ 90fb5f48-5f5c-4fc8-a431-7179ea0a79b2
md"""
#### 6. 行列の特異値分解による低階最良近似
一般の $m{\times}n$ 行列 $A=\begin{bmatrix}\vec{a}_1 & \dotsb & \vec{a}_n \end{bmatrix}$ に対して $A$ の階数 (rank) とよばれる非負整数

$r
:=
\operatorname{dim}\bigl(
\{
A\vec{x}
=
x_1\vec{a}_1+\dotsb+x_n\vec{a}_n 
\ \vert \ 
\vec{x} \in \mathbb{R}^n
\}
\bigr)$ 

が定義される. $A$ の成分がすべて $0$ の場合を除くと $1 \leqq r \leqq \min\{m,n\}$ が成立する.

このとき 

$\sigma_1 \geqq \dotsb \geqq \sigma_r>0,
\quad
\vec{u}_1,\dotsc,\vec{u}_r \in \mathbb{R}^m,
\quad
\vec{v}_1,\dotsc,\vec{v}_r \in \mathbb{R}^n$

が存在して

$A=\sum_{j=1}^r\sigma_j\vec{u}_j\vec{v}_j^T$

と表すことができることが知られている. 

$A_k:=\sum_{j=1}^k\sigma_j\vec{u}_j\vec{v}_j^T, \quad k=1,\dotsc,r$

とすると, $A_k$ の階数は $k$ であり, 階数は $k$ のすべての行列のうち $A$ に最も近い行列, 
すなわち階数は $k$ のすべての行列の中で $A_k$ の最良の近似を与える行列であることが知られている. 
香港城市大学 (City University of Hong Kong) の食堂の経済飯の画像 $384\times512\times3$ を低階最良近似してみる. 
"""

# ╔═╡ 05f8285f-5388-41ec-8c3f-e94a1f91cb2d
md"""
rank = $(@bind r１ Slider(1:rank, show_value=true))
"""

# ╔═╡ 9cc5fb8b-b53f-45d5-a8a9-c2051be8ded9
begin
SVD1=plot(colorview(RGB,W[:,:,:,r１]),
        title="approximation 1-100/384",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
SVD2=plot(colorview(RGB,X),
        title="original RGB image",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
plot(SVD1,SVD2)
end

# ╔═╡ 179fda5c-5ef8-4990-9f37-05b2668a3583


# ╔═╡ 009daab5-c2fd-42bb-addf-995b0f9d7eca
md"""
#### 7. 離散ハール・ウェーブレットとウェーブレット分解
離散ウェーブレットとはある性質をみたす2つの数ベクトルの組でのことである. 
最も単純なのは次のハール・ウェーブレットである. 

$\vec{u}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1
\\
1
\\
0
\\
\vdots
\\
0 
\end{bmatrix}, 
\quad
\vec{v}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
1
\\
-1
\\
0
\\
\vdots
\\
0 
\end{bmatrix}
\in\mathbb{R}^{N},$

これを用いると数ベクトルや行列の隣接する成分の平均をとるという操作を線形代数だけを用いて系統的に定義することができる. 2回の操作で

$\vec{a}_0
=
\begin{bmatrix}
1 \\ 3 \\ 5 \\ 7 
\end{bmatrix}
\mapsto 
\vec{a}_1
=
\begin{bmatrix}
2 \\ 2 \\ 6 \\ 6 
\end{bmatrix}
\mapsto 
\vec{a}_2
=
\begin{bmatrix}
4 \\ 4 \\ 4 \\ 4 
\end{bmatrix}$

$A_0
=
\begin{bmatrix}
0 & 2 & 4 & 6
\\ 
8 & 10 & 12 & 14 
\\ 
16 & 18 & 20 & 22
\\ 
24 & 26 & 28 & 30
\end{bmatrix}
\mapsto 
A_1
=
\begin{bmatrix}
5 & 5 & 9 & 9
\\ 
5 & 5 & 9 & 9 
\\ 
21 & 21 & 25 & 25
\\ 
21 & 21 & 25 & 25
\end{bmatrix}
\mapsto 

A_2
=
\begin{bmatrix}
15 & 15 & 15 & 15
\\ 
15 & 15 & 15 & 15 
\\ 
15 & 15 & 15 & 15
\\ 
15 & 15 & 15 & 15
\end{bmatrix}$

のようにな. 平均をとる操作を $\ell$ 回行って得られる $\vec{a}_\ell$ や $A_\ell$ をレベル $\ell$ の近似部分, 残りの $\vec{a} - \vec{a}_\ell$ や $A - A_\ell$ をレベル $\ell$ の詳細部分という. 

$\vec{a}_0 = \vec{a}_\ell + (\vec{a}_0-\vec{a}_\ell),
\quad
A_0 = A_\ell + (A_0-A_\ell),$

画像データのウェーブレット分解を観察してみよう.

""" 

# ╔═╡ 65f12e5d-d912-4d2d-9148-fe6b5a821f17
begin
	L7=3;
	p7=8;
	q7=8;
    G7=[  0   4   8  12  16  20  24  28;
		  32  36  40  44  48  52  56  60;
		  64  68  72  76  80  84  88  92;
		  96 100 104 108 112 116 120 124;
		 128 132 136 140 144 148 152 156;
		 160 164 168 172 176 180 184 188; 
		 192 196 200 204 208 212 216 220;
		 224 228 232 236 240 244 248 252]/255; 

	XG7=zeros(p7,q7,L7);
	for l=1:L7
		XG7[:,:,l]=dwt(G7, wavelet(WT.haar), l);
	end

	XG7approx=zeros(p7,q7,L7);	
	for l=1:L7
	    XG7approx[1:Int(p7/2^l),1:Int(q7/2^l),l]=XG7[1:Int(p7/2^l),1:Int(q7/2^l),l];
	end
	XG7detail=XG7-XG7approx;

    YG7approx=zeros(p7,q7,L7);
	YG7detail=zeros(p7,q7,L7);
	for l=1:L7
		YG7approx[:,:,l]=idwt(XG7approx[:,:,l], wavelet(WT.haar), l);
		YG7detail[:,:,l]=idwt(XG7detail[:,:,l], wavelet(WT.haar), l);
	end
	
    ZG7approx=zeros(p7,q7,L7+1);
	ZG7detail=zeros(p7,q7,L7+1);
	ZG7approx[:,:,1]=G7;
	for l=2:L7+1
        ZG7approx[:,:,l]=YG7approx[:,:,l-1];
	    ZG7detail[:,:,l]=YG7detail[:,:,l-1];
	end
		
end

# ╔═╡ dfb617db-57e3-4704-8621-fb2b18adad20
md"""
#### $8\times8$ 白黒画像
Level = $(@bind l7 Slider(0:L7, show_value=true))
"""

# ╔═╡ 24277486-9446-4e7a-b4d1-e1fedff8ba1a
begin
Z1=plot(Gray.(G7),
        #title="Original Image",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
Z2=plot(Gray.(ZG7approx[:,:,l7+1]),
        #title="Haar",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
Z3=plot(Gray.(ZG7detail[:,:,l7+1]),
        #title="Daubechies 2",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
plot(Z1,Z2,Z3,layout=(1,3))
end

# ╔═╡ 48b3fdb1-952e-49c4-9879-85dba253d467


# ╔═╡ a554ca47-ed38-432e-8330-5b1b03492856
begin

# Decomposition Filter 
	RXII=zeros(2^L1,2^L1,L1);
	GXII=zeros(2^L1,2^L1,L1);
	BXII=zeros(2^L1,2^L1,L1);
	for l=1:L1
		RXII[:,:,l]=dwt(RII, wavelet(WT.haar), l);
		GXII[:,:,l]=dwt(GII, wavelet(WT.haar), l);
		BXII[:,:,l]=dwt(BII, wavelet(WT.haar), l);
	end

# Splitting 
	RXIIapprox=zeros(2^L1,2^L1,L1);
	GXIIapprox=zeros(2^L1,2^L1,L1);
	BXIIapprox=zeros(2^L1,2^L1,L1);
	for l=1:L1
	    RXIIapprox[1:2^(L1-l),1:2^(L1-l),l]=RXII[1:2^(L1-l),1:2^(L1-l),l];
		GXIIapprox[1:2^(L1-l),1:2^(L1-l),l]=GXII[1:2^(L1-l),1:2^(L1-l),l];
		BXIIapprox[1:2^(L1-l),1:2^(L1-l),l]=BXII[1:2^(L1-l),1:2^(L1-l),l];
	end
	RXIIdetail=RXII-RXIIapprox;
	GXIIdetail=GXII-GXIIapprox;
	BXIIdetail=BXII-BXIIapprox;

# Composition Filter
	RYIIapprox=zeros(2^L1,2^L1,L1);
	RYIIdetail=zeros(2^L1,2^L1,L1);
	GYIIapprox=zeros(2^L1,2^L1,L1);
	GYIIdetail=zeros(2^L1,2^L1,L1);
	BYIIapprox=zeros(2^L1,2^L1,L1);
	BYIIdetail=zeros(2^L1,2^L1,L1);
	for l=1:L1
		RYIIapprox[:,:,l]=idwt(RXIIapprox[:,:,l], wavelet(WT.haar), l);
		RYIIdetail[:,:,l]=idwt(RXIIdetail[:,:,l], wavelet(WT.haar), l);
		GYIIapprox[:,:,l]=idwt(GXIIapprox[:,:,l], wavelet(WT.haar), l);
		GYIIdetail[:,:,l]=idwt(GXIIdetail[:,:,l], wavelet(WT.haar), l);
		BYIIapprox[:,:,l]=idwt(BXIIapprox[:,:,l], wavelet(WT.haar), l);
		BYIIdetail[:,:,l]=idwt(BXIIdetail[:,:,l], wavelet(WT.haar), l);
	end

# RGB
	ZIIapprox=zeros(3,2^L1,2^L1,L1+1);
	ZIIdetail=zeros(3,2^L1,2^L1,L1+1);
	ZIIapprox[1,:,:,1]=RII;
	ZIIapprox[2,:,:,1]=GII;
	ZIIapprox[3,:,:,1]=BII;
	for l=2:L1+1
		ZIIapprox[1,:,:,l]=RYIIapprox[:,:,l-1];
		ZIIapprox[2,:,:,l]=GYIIapprox[:,:,l-1];
		ZIIapprox[3,:,:,l]=BYIIapprox[:,:,l-1];
		ZIIdetail[1,:,:,l]=RYIIdetail[:,:,l-1];
		ZIIdetail[2,:,:,l]=GYIIdetail[:,:,l-1];
		ZIIdetail[3,:,:,l]=BYIIdetail[:,:,l-1];
	end

end

# ╔═╡ 57406b4f-ceed-4e0b-96ad-15db166defae
md"""
#### 乱数で生成した $16\times16$ RGB画像
Level = $(@bind l8 Slider(0:L1, show_value=true))
"""

# ╔═╡ cb3d9fe0-1245-4918-b4d5-9caad49fcab5
begin
Q0=plot(colorview(RGB,ZIIapprox[:,:,:,1]),
        #title="Original Image",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
Q1=plot(colorview(RGB,ZIIapprox[:,:,:,l8+1]),
        #title="Approximation",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false);
Q2=plot(colorview(RGB,ZIIdetail[:,:,:,l8+1]),
        #title="Detail",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
plot(Q0,Q1,Q2,layout=(1,3))
end

# ╔═╡ 4c4e797e-0b04-476b-89f8-0a0b118c80ab
begin

# Decomposition Filter 
	XR=zeros(p,q,L);
	XG=zeros(p,q,L);
	XB=zeros(p,q,L);
	for l=1:L
		XR[:,:,l]=dwt(R, wavelet(WT.haar), l);
		XG[:,:,l]=dwt(G, wavelet(WT.haar), l);
		XB[:,:,l]=dwt(B, wavelet(WT.haar), l);
	end

# Splitting 
	XRapprox=zeros(p,q,L);
	XGapprox=zeros(p,q,L);
	XBapprox=zeros(p,q,L);	
	for l=1:L
	XRapprox[1:Int(p/2^l),1:Int(q/2^l),l]=XR[1:Int(p/2^l),1:Int(q/2^l),l];
    XGapprox[1:Int(p/2^l),1:Int(q/2^l),l]=XG[1:Int(p/2^l),1:Int(q/2^l),l];
	XBapprox[1:Int(p/2^l),1:Int(q/2^l),l]=XB[1:Int(p/2^l),1:Int(q/2^l),l];
	end
	XRdetail=XR-XRapprox;
	XGdetail=XG-XGapprox;
	XBdetail=XB-XBapprox;

# Composition Filter
	YRapprox=zeros(p,q,L);
	YGapprox=zeros(p,q,L);
	YBapprox=zeros(p,q,L);
	YRdetail=zeros(p,q,L);
	YGdetail=zeros(p,q,L);
	YBdetail=zeros(p,q,L);
	for l=1:L
		YRapprox[:,:,l]=idwt(XRapprox[:,:,l], wavelet(WT.haar), l);
		YGapprox[:,:,l]=idwt(XGapprox[:,:,l], wavelet(WT.haar), l);
		YBapprox[:,:,l]=idwt(XBapprox[:,:,l], wavelet(WT.haar), l);
		YRdetail[:,:,l]=idwt(XRdetail[:,:,l], wavelet(WT.haar), l);
		YGdetail[:,:,l]=idwt(XGdetail[:,:,l], wavelet(WT.haar), l);
		YBdetail[:,:,l]=idwt(XBdetail[:,:,l], wavelet(WT.haar), l);
	end
	
# RGB
	Wapprox=zeros(3,p,q,L+1);
	Wdetail=zeros(3,p,q,L+1);
	Wapprox[:,:,:,1]=A;
	for l=2:L+1
		Wapprox[1,:,:,l]=YRapprox[:,:,l-1];
		Wapprox[2,:,:,l]=YGapprox[:,:,l-1];
		Wapprox[3,:,:,l]=YBapprox[:,:,l-1];
		Wdetail[1,:,:,l]=YRdetail[:,:,l-1];
		Wdetail[2,:,:,l]=YGdetail[:,:,l-1];
		Wdetail[3,:,:,l]=YBdetail[:,:,l-1];
	end
end


# ╔═╡ f5961a2c-8ebc-4fb8-bb25-108372588c38
md"""
#### $384\times512$ RGB画像
##### 香港城市大学の食堂の経済飯 $384=2^7\times3$, $512=2^9$
Level = $(@bind l Slider(0:L, show_value=true))
"""

# ╔═╡ 84911507-1548-445a-942c-b0ca93cab0b6
begin
P1=plot(colorview(RGB,A[:,:,:]),
        #title="Approximation",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
P2=plot(colorview(RGB,Wapprox[:,:,:,l+1]),
        #title="Approximation",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
P3=plot(colorview(RGB,Wdetail[:,:,:,l+1]),
        #title="Detail",
        xaxis=false, 
        xticks=false, 
        yaxis=false, 
        yticks=false, 
        grid=false);
plot(P1,P2,P3,size=(1536,384),layout=(1,3))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Formatting = "59287772-0a20-5a39-b81b-1366585eb4c0"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
ImageMagick = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
ImageTransformations = "02fcd773-0e25-5acc-982a-7f6622650795"
ImageView = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
ImplicitPlots = "55ecb840-b828-11e9-1645-43f4a9f9ace7"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LowRankApprox = "898213cb-b102-5a47-900c-97e73b919f73"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Primes = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TestImages = "5e47fb64-e119-507b-a336-dd2b206d9990"
Wavelets = "29a6e085-ba6d-5f35-a997-948ac2efa89a"

[compat]
Colors = "~0.13.0"
DSP = "~0.8.3"
DataFrames = "~1.7.0"
Formatting = "~0.4.3"
ForwardDiff = "~0.10.36"
ImageFiltering = "~0.7.9"
ImageMagick = "~1.4.1"
ImageTransformations = "~0.10.1"
ImageView = "~0.12.6"
Images = "~0.26.2"
ImplicitPlots = "~0.2.3"
LaTeXStrings = "~1.4.0"
LowRankApprox = "~0.5.5"
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
Primes = "~0.5.6"
Symbolics = "~6.38.0"
TestImages = "~1.9.0"
Wavelets = "~0.10.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "b5b97ce57e86b03028dff3a3907ca3b3678e8107"

[[deps.ADTypes]]
git-tree-sha1 = "e2478490447631aedba0823d4d7a80b2cc8cdb32"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.14.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "0c8a70a69036d8f3a2426d768d30144547cd73c0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.3"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "6d8b535fd38293bc54b88455465a1386f8ac1c3c"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.119"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "9a3ae38b460449cc9e7dd0cfb059c76028724627"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.1"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8e2d86e06ceb4580110d9e716be26658effc5bfd"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "da121cbdc95b065da07fbb93638367737969693f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.8+0"

[[deps.GTK4_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "Graphene_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Libepoxy_jll", "Pango_jll", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libX11_jll", "Xorg_libXcursor_jll", "Xorg_libXdamage_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "Xorg_libXrender_jll", "gdk_pixbuf_jll", "iso_codes_jll", "xkbcommon_jll"]
git-tree-sha1 = "acca2e4afa083ca045b9ccef613d760d9e6d30bf"
uuid = "6ebb71f1-8434-552f-b6b1-dc18babcca63"
version = "4.12.4+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphene_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "61850a17f562453e3485a489c9c8cccb3abcab93"
uuid = "75302f13-0b7e-5bab-a6d1-23fa92e4c2ea"
version = "1.10.6+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "3169fd3440a02f35e549728b0890904cfd4ae58a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Gtk4]]
deps = ["BitFlags", "CEnum", "Cairo", "Cairo_jll", "ColorTypes", "FixedPointNumbers", "GTK4_jll", "Glib_jll", "Graphene_jll", "Graphics", "JLLWrappers", "Libdl", "Librsvg_jll", "Pango_jll", "Preferences", "Reexport", "Scratch", "Xorg_xkeyboard_config_jll", "adwaita_icon_theme_jll", "gdk_pixbuf_jll", "hicolor_icon_theme_jll", "libpng_jll"]
git-tree-sha1 = "786eb3183b7ca7bf57b24a0e07cbb0f0de9bc876"
uuid = "9db2cae5-386f-4011-9d63-a5602296539b"
version = "0.6.8"

[[deps.GtkObservables]]
deps = ["Cairo", "Colors", "Dates", "FixedPointNumbers", "Graphics", "Gtk4", "IntervalSets", "LinearAlgebra", "Observables", "PrecompileTools", "Reexport", "RoundingIntegers"]
git-tree-sha1 = "00111abceb335becbbc4d38f26228671f5dd6ecc"
uuid = "8710efd8-4ad6-11eb-33ea-2d5ceb25a41c"
version = "2.1.4"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "f93655dc73d7a0b4a368e3c0bce296ae035ad76e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.16"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "33cb509839cc4011beb45bde2316e64344b0f92b"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.9"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "8582eca423c1c64aac78a607308ba0313eeaed56"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.4.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fa01c98985be12e5d75301c4527fff2c46fa3e0e"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "7.1.1+1"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "cffa21df12f00ca1a365eb8ed107614b40e8c6da"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.6"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "3db3bb9f7014e86f13692581fa2feb6460bdee7e"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.4"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "e0884bdf01bbbb111aea77c348368a86fb4b5ab6"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.1"

[[deps.ImageView]]
deps = ["AxisArrays", "Cairo", "Compat", "Graphics", "Gtk4", "GtkObservables", "ImageBase", "ImageCore", "ImageMetadata", "MultiChannelColors", "PrecompileTools", "RoundingIntegers", "StatsBase"]
git-tree-sha1 = "e535c709a4fb6f0ae65026fa9bfac624abec016f"
uuid = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
version = "0.12.6"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "a49b96fd4a8d1a9a718dfd9cde34c154fc84fcd5"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.ImplicitPlots]]
deps = ["Contour", "MultivariatePolynomials", "RecipesBase", "Requires", "StaticArrays", "StaticPolynomials"]
git-tree-sha1 = "baaa32fec0346ccf55b61972858f33809d8f9694"
uuid = "55ecb840-b828-11e9-1645-43f4a9f9ace7"
version = "0.2.3"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "b842cbff3f44804a84fda409745cc8f04c029a20"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "8e071648610caa2d3a5351aba03a936a0c37ec61"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.13"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libepoxy_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "7a0158b71f8be5c771e7a273183b2d0ac35278c5"
uuid = "42c93a91-0102-5b3f-8f9d-e41de60ac950"
version = "1.5.10+0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "e5afce7eaf5b5ca0d444bcb4dc4fd78c54cbbac0"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.172"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LowRankApprox]]
deps = ["FFTW", "LinearAlgebra", "LowRankMatrices", "Nullables", "Random", "SparseArrays"]
git-tree-sha1 = "031af63ba945e23424815014ba0e59c28f5aed32"
uuid = "898213cb-b102-5a47-900c-97e73b919f73"
version = "0.5.5"

[[deps.LowRankMatrices]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7c8664b2f3d5c3d9b77605c03d53b18813e79b0f"
uuid = "e65ccdef-c354-471a-8090-89bec1c20ec3"
version = "1.0.1"
weakdeps = ["FillArrays"]

    [deps.LowRankMatrices.extensions]
    LowRankMatricesFillArraysExt = "FillArrays"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "e9650bea7f91c3397eb9ae6377343963a22bf5b8"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.8.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "453de0fc2be3d11b9b93ca4d0fddd91196dcf1ed"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.5"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultiChannelColors]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "Compat", "FixedPointNumbers", "LinearAlgebra", "Reexport", "Requires"]
git-tree-sha1 = "c4dce3e565ee81dccd3b588a7ea08cfc67778339"
uuid = "d4071afc-4203-49ee-90bc-13ebeb18d604"
version = "0.1.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "8a3271d8309285f4db73b4f662b1b290c715e85e"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.21"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Nullables]]
git-tree-sha1 = "8f87854cc8f3685a60689d8edecaa29d2251979b"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "1.0.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ad31332567b189f508a3ea8957a2640b1147ab00"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "0e1340b5d98971513bddaa6bbed470670cebbbfe"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.34"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "44f6c1f38f77cafef9450ff93946c53bd9ca16ff"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "f202a1ca4f6e165238d8175df63a7e26a51e04dc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.7"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "555c272d20fc80a2658587fb9bbda60067b93b7c"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.19"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "112c876cee36a5784df19098b55db2b238afc36a"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.31.2"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RoundingIntegers]]
git-tree-sha1 = "99acd97f396ea71a5be06ba6de5c9defe188a778"
uuid = "d5f540fe-1c90-5db3-b776-2e2f362d9394"
version = "1.1.0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "1f7cf417da3771b98f0e3f32ce0bb813e9fe91fa"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.85.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "1c4b7f6c3e14e6de0af66e66b86d525cae10ecb4"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.13"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.StaticPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "StaticArrays"]
git-tree-sha1 = "0b4ec86a5ba2269c51897381dd0e7a222650f447"
uuid = "62e018b1-6e46-5407-a5a7-97d4fbcae734"
version = "1.3.7"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "35b09e80be285516e52c9054792c884b9216ae3c"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.4.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringDistances]]
deps = ["Distances", "StatsAPI"]
git-tree-sha1 = "5b2ca70b099f91e54d98064d5caf5cc9b541ad06"
uuid = "88034a9c-02f8-509d-84a9-84ec65e18404"
version = "0.11.3"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "7530e17b6ac652b009966f8ad53371a4ffd273f2"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.39"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "2c9879cd67d1bb2f2989669e5849639bb4d3c792"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.26.1"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "e46dbf646bc3944c22a37745361c2e0a94f81d66"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.38.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "d155450e6dff2a8bc2fcb81dcb194bd98b0aeb46"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.2"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestImages]]
deps = ["AxisArrays", "ColorTypes", "FileIO", "ImageIO", "ImageMagick", "OffsetArrays", "Pkg", "StringDistances"]
git-tree-sha1 = "fc32a2c7972e2829f34cf7ef10bbcb11c9b0a54c"
uuid = "5e47fb64-e119-507b-a336-dd2b206d9990"
version = "1.9.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "18ad3613e129312fe67789a71720c3747e598a61"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.3"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "f21231b166166bebc73b99cea236071eb047525b"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.3"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f57facfd1be61c42321765d3551b3df50f7e09f6"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.28"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "4ab62a49f1d8d9548a1c8d1a75e5f55cf196f64e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.71"

[[deps.Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d0ec97a100abbe47a5e9a02361841da49cce6029"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.10.1"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5db3e9d307d32baba7067b13fc7b5aa6edd4a19a"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.36.0+0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "807c226eaf3651e7b2c468f687ac788291f9a89b"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.3+0"

[[deps.Xorg_libXdamage_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll"]
git-tree-sha1 = "534ed7d299469f3438b2c136d7beb0b50da88ce0"
uuid = "0aeada51-83db-5f97-b67e-184615cfc6f6"
version = "1.1.6+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "984b313b049c89739075b8e2a94407076de17449"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.2+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.adwaita_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "hicolor_icon_theme_jll"]
git-tree-sha1 = "28401767f30e5743ef5e3b0be71417bc911d3952"
uuid = "b437f822-2cd6-5e08-a15c-8bac984d38ee"
version = "43.0.1+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "e9190f9fb03f9c3b15b9fb0c380b0d57a3c8ea39"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.8+0"

[[deps.hicolor_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b458a6f6fc2b1a8ca74ed63852e4eaf43fb9f5ea"
uuid = "059c91fe-1bad-52ad-bddd-f7b78713c282"
version = "0.17.0+3"

[[deps.iso_codes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4d295b7797afbe24ad5d0452f673b901a81c43c3"
uuid = "bf975903-5238-5d20-8243-bc370bc1e7e5"
version = "4.17.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "ccbb625a89ec6195856a50aa2b668a5c08712c94"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.4.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "63406453ed9b33a0df95d570816d5366c92b7809"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+2"
"""

# ╔═╡ Cell order:
# ╟─5ffc5ec0-f91f-11ec-0552-37ef5f25102d
# ╟─f3034480-e17d-4a7b-bdf3-cfb1d55d7bb9
# ╟─0a085483-72c0-4ec5-a677-9cf09f05cc54
# ╟─2aa2a4f4-65d4-43e4-b8fb-13bcca1b84b9
# ╟─7d7e3384-cce9-4a97-bd78-604932418d56
# ╟─38a19ae1-07b9-414c-999c-c74d2b8df7d5
# ╟─c2b971bd-5dc2-4f84-80d2-dd5307a96961
# ╟─763ca8d2-4975-4599-9aa2-247d2775c0cc
# ╟─84e0ec88-3824-4ad8-a38c-c01d9cacc7d7
# ╟─5b8d6598-f602-4408-b948-84cec2bfc52a
# ╟─05c00c05-c990-4632-ba43-4f6e14fad19f
# ╟─4c26a3e2-7ad7-4ca5-a5bb-1cf747f1a11a
# ╟─769adba4-34e2-4e63-9099-2df3f6dbb529
# ╟─c0233b8e-9db1-4a3a-8080-9f4bf7f3804b
# ╟─76cb2ec2-fc76-4ad5-8d57-777395755f7f
# ╟─8ed1312b-f5a3-4220-a865-9a12ed09a691
# ╟─21c924bb-23de-42ac-805b-a19744a8e105
# ╟─cfb114a5-e5c4-429d-bf6e-9acf7ea4b31d
# ╟─3c696822-5b77-4512-83f0-c524371cde5d
# ╟─728a353b-e2fe-460d-8dcc-0f4cadbbd033
# ╟─2db362e9-e32a-4ba2-8810-3bdef318af41
# ╟─6d2d8054-56e3-4499-a8ea-cfcf291361e8
# ╟─b61b5be5-2b6e-4450-a3f1-6abe61bf9974
# ╟─2dc43448-3cf4-475f-8b61-d16c3594d9e3
# ╟─06191b8b-ddcc-4a2b-8228-0680aa8ad7c5
# ╟─b5e407b0-f3b2-4c9c-9274-9f4a69660110
# ╟─f64936c9-a313-4c4b-a0b1-5e0df1b9e508
# ╟─dbc102fd-f90e-4083-9155-173a0fd090a5
# ╟─21fda1bd-15c8-40e9-bb4f-ebbbd6c102ef
# ╟─7bb3421a-89ab-4566-aafc-c07350243d9c
# ╟─df9b6dd3-9e44-4479-8049-233c774d345d
# ╟─0ce69d78-5b11-4012-91e2-cbb9bde9e22f
# ╟─4d7243c8-126d-40ac-a579-e53406a102d9
# ╟─6efa7afe-b600-4ce4-81cb-43aa7348d8c2
# ╟─4edf5113-e3a3-4f00-9fbc-1aebe02dfed1
# ╟─baa5dd97-3152-469f-a45e-d427c786dcfd
# ╟─407adc66-f13c-4de1-8ad8-f59bf6b66063
# ╟─71e8ddf6-8c83-46ad-9222-3f1a9f57f2bf
# ╟─352fa6e4-1fc2-48ba-a62f-0e60ca0a6304
# ╟─76fa77f0-d82e-40ac-b826-50b8c781d030
# ╟─b742bed1-9ae5-44f7-8b32-6609cdc3741b
# ╟─6facfa26-5ae1-484f-a4f0-e8846ee8aec4
# ╟─50841743-af5c-4cea-9e61-6514ee36ec72
# ╟─bd09658a-e155-42ef-bebe-5088f55f177f
# ╟─7fbec145-da40-4636-875d-c5124c6da60b
# ╟─d9b4f11f-3687-4872-8237-bcfb396e9d88
# ╟─14b5204e-f53d-4f5c-aaa1-070dde6eed35
# ╟─b469f032-77fa-48a3-976d-79de94ddea4a
# ╟─5b538799-ea23-46cf-9f40-34e902cfc91d
# ╟─2196331a-3546-47f8-a423-5c49441e2bbb
# ╟─4843d19a-ef5d-47d5-8c23-8dedeb069641
# ╟─feaf5998-8014-4ff6-b600-70af413ebbba
# ╟─d2d1002e-12d0-4470-8681-731003de6c6d
# ╟─431d018d-c008-47f3-9d30-7e9da54d579a
# ╟─8d9148e7-4ce1-4da7-9eb3-43df098d5ad5
# ╟─ab1c3d69-0499-4fad-b2a8-5360ed7339e1
# ╟─0e933e40-e613-4eda-bc8c-666a7d208c02
# ╟─b6af8639-aaf1-4148-9ddc-410d2a7b0495
# ╟─c15d7baf-54d5-49da-b144-3ce7255cc2f1
# ╟─3f17e529-3ca4-4302-a54b-30fe096893ac
# ╟─157db428-8868-403b-9e3b-705f3df7612a
# ╟─f95d5371-c27e-459f-9ca8-b36c4aa67177
# ╟─fbad5a13-c986-4d2a-936d-51757d7eac93
# ╟─1299f123-8d64-423f-bc64-6410a7392cc2
# ╟─1c624c90-99ed-414d-a6e4-b2f6b16101c8
# ╟─0969f38b-f1f6-4b81-b358-c87ab486beff
# ╟─3e57202c-4c53-4282-a477-410a990bb318
# ╟─659ab5e5-5b78-4496-a3d5-5f4b22396178
# ╟─ca6d79ef-d5a4-4a55-b7f9-83d8682ddd1b
# ╟─b06741a8-5ca6-4e15-b085-e3bbcddca93d
# ╟─87dd79c0-c0e9-4021-879a-c18cfa6e9252
# ╟─f68a62e7-f442-49d1-bd1d-b7268d55f037
# ╟─3e301f5d-cdb8-43c2-8562-a83416c270cb
# ╟─173a1432-4ca9-4ce1-9de2-2ef441b1810d
# ╟─3cd15536-0761-44db-a8d1-96096d24e7c5
# ╟─213f0010-d7b9-4584-8484-868c5cd9316d
# ╟─f7a1f0f0-2cd1-412b-98e2-9b3bb71b5c43
# ╟─295e359a-ae03-4d15-b6f9-ee7c98067eb9
# ╟─3315d5ae-444e-4958-a1d7-086eeaaad522
# ╟─7b45973d-bfce-4b9b-8d2c-e35b5ef421a5
# ╟─2a819b34-8fdc-47b9-9397-72e8f791a902
# ╟─6450807b-e3e9-43e6-bd64-50c385c077ce
# ╟─fe18cdd8-c5d6-443f-b99d-66514a9393b8
# ╟─1de04ede-fa17-4497-a85d-2205df043d04
# ╟─0f18af6e-c2ba-43df-8299-2c42821c8c17
# ╟─de50b16d-e218-4221-ba1a-a4539dc33973
# ╟─02d137fc-4aac-4da4-8850-04a6529b7042
# ╟─9f9c9544-7454-4277-84a4-7a2b46fdea2f
# ╟─8c938021-1a02-44f6-93b0-a474f89a60d4
# ╟─8c14a644-b7f3-45eb-bc0c-bc9fc28d65bd
# ╟─4d3d1122-2424-478f-894d-f7180e59fca2
# ╟─8ceef78f-5058-4e6a-b5ce-ea5b22f6a7fb
# ╟─83b7d146-00a0-47b7-8988-288721afad10
# ╟─90fb5f48-5f5c-4fc8-a431-7179ea0a79b2
# ╟─05f8285f-5388-41ec-8c3f-e94a1f91cb2d
# ╟─9cc5fb8b-b53f-45d5-a8a9-c2051be8ded9
# ╟─179fda5c-5ef8-4990-9f37-05b2668a3583
# ╟─009daab5-c2fd-42bb-addf-995b0f9d7eca
# ╟─65f12e5d-d912-4d2d-9148-fe6b5a821f17
# ╟─dfb617db-57e3-4704-8621-fb2b18adad20
# ╟─24277486-9446-4e7a-b4d1-e1fedff8ba1a
# ╟─48b3fdb1-952e-49c4-9879-85dba253d467
# ╟─a554ca47-ed38-432e-8330-5b1b03492856
# ╟─57406b4f-ceed-4e0b-96ad-15db166defae
# ╟─cb3d9fe0-1245-4918-b4d5-9caad49fcab5
# ╟─4c4e797e-0b04-476b-89f8-0a0b118c80ab
# ╟─f5961a2c-8ebc-4fb8-bb25-108372588c38
# ╟─84911507-1548-445a-942c-b0ca93cab0b6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
