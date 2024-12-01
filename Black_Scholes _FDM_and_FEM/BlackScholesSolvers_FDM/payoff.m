function result = payoff(stockPrice, strikePrice, isCallOption)
    if isCallOption
        result = max(stockPrice - strikePrice, 0);
    else
        result = max(strikePrice - stockPrice, 0);
    end
    end
    