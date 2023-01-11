function [] = print_msg(msg_priority,msg,ctrl)
	% print mesagge into stdout according to priority
	if (msg_priority <= ctrl.verbose)
		fprintf('%s\n',msg)
	end
