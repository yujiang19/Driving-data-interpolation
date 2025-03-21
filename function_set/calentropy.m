function entropy =  calentropy(pdf_values, x_values, f_max)
indices = pdf_values > 0;
pdf_values = pdf_values(indices);
x_values = x_values(indices);
p = pdf_values./sum(pdf_values);
dx = abs(x_values(2) - x_values(1)); 
entropy = -sum(p .* log2(p)).*dx/f_max;  
end