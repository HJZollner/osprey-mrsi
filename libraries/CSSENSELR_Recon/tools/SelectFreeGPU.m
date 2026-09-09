function thisgpu = SelectFreeGPU(MemoryNeeded)
	Ngpu=gpuDeviceCount();%"available");
	if Ngpu>0
		for gpu=1:Ngpu
			gpuDevice([]);
			D=gpuDevice(gpu);
			GPUMemory(gpu)=D.AvailableMemory;
		end
		[GPUAvailableMemory,thisgpu]=max(GPUMemory);	
		gpuDevice([]); gpuDevice(thisgpu);	
		
		try
		    nnet.internal.cnngpu.reluForward(1);
		catch ME
		end
		
		if GPUAvailableMemory<MemoryNeeded
			warning('No GPU with required free memory could be found! Running on CPU!');	
			thisgpu=-1;
		end	
	else
		warning('No GPU available to Matlab on this system. Running on CPU!');
		thisgpu=-1;
	end

end
