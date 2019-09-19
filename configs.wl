cosine = Cos[2 * Pi / 7];
ch2phi = (8/3 * Cos[Pi / 7]^2) - 1;
rsquare = (1 - cosine)/(ch2phi - 1);
x0square = rsquare + 1;
r = Sqrt[rsquare];
x0 = FullSimplify[Sqrt[x0square]];

configs = <|
	"h_7_3" -> {
		{
			{0, 0, 1},
			{0, -Sin[Pi/7], Cos[Pi/7]},
			{Sqrt[3]/2, Sqrt[3/4 + rsquare], 0}
		}, 
		{x0, Cos[Pi/7], Sin[Pi/7]},
		Red
	},

	"h_3_7" -> {
		{
			{Sqrt[3]/2, Sqrt[3/4 + rsquare], 0},
			{0, -Sin[Pi/7], Cos[Pi/7]},
			{0, 0, 1}
		},
		{r, 0, 0},
		Cyan
	},

	"h_72_7" -> {
		{
			{0, -Sin[2 Pi/7], Cos[2 Pi/7]},
			{0, 0, 1},
			{Sqrt[3]/2, Sqrt[3/4 + rsquare], 0}
		}, 
		{Root[-27 - 756*#1^2 - 2688*#1^4 + 448*#1^6 & , 2, 0], Sqrt[3/(2 - 2*Sin[(3*Pi)/14])], 0},
		Blue
	},

	"h_7_72" -> {
		{
			{0, -Sin[Pi/7], Cos[Pi/7]},
			{0, 0, 1},
			{Sqrt[3]/2, Sqrt[3/4 + rsquare], 0}
		}, 
		{Root[-27 - 756*#1^2 - 2688*#1^4 + 448*#1^6 & , 2, 0], Sqrt[3/(2 - 2*Sin[(3*Pi)/14])], 0},
		Yellow
	},

	"null" -> {}
|>
