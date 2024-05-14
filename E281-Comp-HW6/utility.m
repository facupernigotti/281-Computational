function u = utility(c, gamma)
    if gamma == 1
        u = log(c);
    else
        u = c.^(1 - gamma) / (1 - gamma);
    end
end
