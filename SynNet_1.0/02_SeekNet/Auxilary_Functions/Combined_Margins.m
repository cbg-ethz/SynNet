function MixedMargin = Combined_Margins(MarginStats)
global Consts
MixedMargin = Consts.Learning_Worst_Margin_Weight * MarginStats.MarginW + Consts.Learning_Avrage_Margin_Weight * MarginStats.MarginA;
end


