function varargout = mtspof_ga(varargin)
    %
    % 初始化默认配置
    %
    defaultConfig.xy          = 10*rand(40,2); % 默认的城市位置，40个随机2D点
    defaultConfig.dmat        = [];            % 城市间的距离矩阵，默认为空
    defaultConfig.nSalesmen   = 5;             % 销售员数量，默认为5
    defaultConfig.minTour     = 1;             % 每个销售员的最小巡回长度，默认为1
    defaultConfig.popSize     = 80;            % 种群大小，默认为80
    defaultConfig.numIter     = 5e3;           % 迭代次数，默认为5000
    defaultConfig.showProg    = true;          % 是否显示进度，默认为真
    defaultConfig.showStatus  = true;          % 是否显示状态，默认为真
    defaultConfig.showResult  = true;          % 是否显示结果，默认为真
    defaultConfig.showWaitbar = false;         % 是否显示等待条，默认为假

    %
    % 解释用户配置输入
    %
    if ~nargin
        userConfig = struct(); % 如果没有输入参数，则创建空的结构体
    elseif isstruct(varargin{1})
        userConfig = varargin{1}; % 如果第一个参数是结构体，则直接使用
    else
        try
            userConfig = struct(varargin{:}); % 尝试将参数/值对转换为结构体
        catch
            error('??? Expected inputs are either a structure or parameter/value pairs');
        end
    end

    %
    % 用用户输入覆盖默认配置
    %
    configStruct = get_config(defaultConfig,userConfig);

    %
    % 提取配置
    %
    xy          = configStruct.xy;
    dmat        = configStruct.dmat;
    nSalesmen   = configStruct.nSalesmen;
    minTour     = configStruct.minTour;
    popSize     = configStruct.popSize;
    numIter     = configStruct.numIter;
    showProg    = configStruct.showProg;
    showStatus  = configStruct.showStatus;
    showResult  = configStruct.showResult;
    showWaitbar = configStruct.showWaitbar;
    if isempty(dmat) % 如果距离矩阵为空，则根据城市位置计算
        nPoints = size(xy,1);
        a = meshgrid(1:nPoints);
        dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
    end

    %
    % 验证输入
    %
    [N,dims] = size(xy); % N是城市数量，dims是维度
    [nr,nc] = size(dmat); % nr和nc是距离矩阵的行数和列数
    if (N ~= nr) || (N ~= nc)
        error('??? Invalid XY or DMAT inputs'); % 如果输入无效，抛出错误
    end
    n = N - 2; % 去除起点和终点后的城市数量

    %
    % 健全性检查
    %
    nSalesmen   = max(1,min(n,round(real(nSalesmen(1)))));
    minTour     = max(1,min(floor(n/nSalesmen),round(real(minTour(1)))));
    popSize     = max(8,8*ceil(popSize(1)/8)); % 种群大小至少为8，且为8的倍数
    numIter     = max(1,round(real(numIter(1)))); % 迭代次数至少为1
    showProg    = logical(showProg(1)); % 转换为逻辑值
    showStatus  = logical(showStatus(1));
    showResult  = logical(showResult(1));
    showWaitbar = logical(showWaitbar(1));

    %
    % 初始化路线分割点选择
    %
    nBreaks = nSalesmen-1; % 分割点数量，即销售员数量减1
    dof = n - minTour*nSalesmen;          % 自由度，即城市数量减去最小巡回长度乘以销售员数量
    addto = ones(1,dof+1);
    for k = 2:nBreaks
        addto = cumsum(addto); % 累加，用于计算分割点的概率分布
    end
    cumProb = cumsum(addto)/sum(addto); % 累积概率

    %
    % 初始化种群
    %
    popRoute = zeros(popSize,n);         % 种群的路线，初始化为零矩阵
    popBreak = zeros(popSize,nBreaks);   % 种群的分割点，初始化为零矩阵
    popRoute(1,:) = (1:n) + 1; % 第一个个体的路线，简单地设置为1到n
    popBreak(1,:) = rand_breaks(); % 第一个个体的分割点，随机生成
    for k = 2:popSize % 为剩下的个体生成随机路线和分割点
        popRoute(k,:) = randperm(n) + 1;
        popBreak(k,:) = rand_breaks();
    end

    %
    % 如果可用，用之前的结果种子算法
    %
    if all(isfield(userConfig,{'optRoute','optBreak'}))
        optRoute = userConfig.optRoute; % 最优路线
        optBreak = userConfig.optBreak; % 最优分割点
        isValidRoute = isequal(popRoute(1,:),sort(optRoute)); % 验证路线是否有效
        isValidBreak = all(optBreak > 0) && all(optBreak <= n) && ...
            (length(optBreak) == nBreaks) && ~any(mod(optBreak,1)); % 验证分割点是否有效
        if isValidRoute && isValidBreak
            popRoute(1,:) = optRoute; % 使用最优路线和分割点初始化第一个个体
            popBreak(1,:) = optBreak;
        end
    end

    %
    % 为绘制的路线选择颜色
    %
    pclr = ~get(0,'DefaultAxesColor'); % 背景颜色的反色
    clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0]; % 默认颜色集
    if (nSalesmen > 5)
        clr = hsv(nSalesmen); % 如果销售员数量超过5，使用HSV颜色映射
    end

    %
    % 运行遗传算法
    %
    row = zeros(popSize,n+nSalesmen); % 初始化用于计算总距离的临时变量
    col = zeros(popSize,n+nSalesmen);
    isValid = false(1,n+nSalesmen); % 标记有效的路线和分割点
    globalMin = Inf; % 全局最小距离
    distHistory = NaN(1,numIter); % 距离历史记录
    tmpPopRoute = zeros(8,n); % 临时种群的路线
    tmpPopBreak = zeros(8,nBreaks); % 临时种群的分割点
    newPopRoute = zeros(popSize,n); % 新一代的路线
    newPopBreak = zeros(popSize,nBreaks); % 新一代的分割点
    [isClosed,isStopped,isCancelled] = deal(false); % 标记变量，用于控制算法的停止和取消
    if showProg % 如果显示进度
        hFig = figure('Name','MTSPOF_GA | Current Best Solution', ...
            'Numbertitle','off','CloseRequestFcn',@close_request); % 创建图形窗口
        hAx = gca; % 获取当前坐标轴句柄
        if showStatus % 如果显示状态
            [hStatus,isCancelled] = figstatus(0,numIter,[],hFig); % 创建状态栏
        end
    end
    if showWaitbar % 如果显示等待条
        hWait = waitbar(0,'Searching for near-optimal solution ...', ...
            'CreateCancelBtn',@cancel_search); % 创建等待条
    end
    isRunning = true; % 标记算法正在运行
    for iter = 1:numIter % 迭代开始
        
        % 此处省略评估解决方案、选择最佳解、修改种群等详细步骤的注释
        
    end
    if showProg && showStatus % 如果显示进度和状态
        figstatus(numIter,numIter,hStatus,hFig); % 更新状态到完成
    end
    if showWaitbar % 如果显示等待条
        delete(hWait); % 关闭等待条
    end
    isRunning = false; % 标记算法停止运行
    if isClosed % 如果窗口已关闭
        delete(hFig); % 删除图形窗口
    end

    % 此处省略返回输出、生成随机分割点、取消搜索和关闭窗口功能的详细注释

end
