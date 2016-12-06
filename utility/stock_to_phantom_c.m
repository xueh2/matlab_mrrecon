function v = stock_to_phantom_c(c_stock, v_phantom, v_stock)
% from the concentration of stock, compute concentration of phantom

v_total = v_phantom+v_stock;
n = c_stock*v_stock;
v = n/v_total;