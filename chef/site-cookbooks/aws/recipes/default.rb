#
# Cookbook Name:: aws
# Recipe:: default
#
# Copyright 2016, GigaScience
#

include_recipe 'user'
include_recipe 'iptables'
#include_recipe 'fail2ban'
include_recipe 'selinux'

# Locates CUDDEL files in /vagrant directory
vagrant_dir = node[:cuddel][:vagrant_dir]

############################
#### Configure iptables ####
############################

service "iptables" do
    action :start
end

iptables_rule 'prefix'
iptables_rule 'http'
iptables_rule 'ssh'
iptables_rule 'postfix'

############################
#### Configure fail2ban ####
############################

# Protect against DDoS attacks
#file "/etc/fail2ban/jail.local" do
#    content <<-EOS
#[ssh-ddos]

#enabled  = true
#port     = ssh
#filter   = sshd-ddos
#logpath  = /var/log/secure
#maxretry = 6
#    EOS
#    owner "root"
#    group "root"
#    mode 0644
#    notifies :restart, "service[fail2ban]"
#end

###########################
#### Configure SELinux ####
###########################

selinux_state 'permissive' do
    action :permissive
end

# Add SELinux policies for GigaBlog
# bash 'gigadb-admin group permissions' do
#     code <<-EOH
#         semanage fcontext -a -t httpd_sys_rw_content_t "/vagrant/log(/.*)?"
#         restorecon -Rv "/vagrant/log"
#
#         semanage fcontext -a -t httpd_sys_rw_content_t "/vagrant/theme(/.*)?"
#                 restorecon -Rv "/vagrant/theme"
#
#         setsebool -P httpd_can_network_connect 1
#         setsebool -P httpd_can_network_connect_db 1
#     EOH
# end

#selinux_state 'enforcing' do
#    action :enforcing
#end

##############################
#### User and group admin ####
##############################

# Create user accounts
user_account node[:cuddel][:user1] do
    comment   node[:cuddel][:user1_name]
    ssh_keys  node[:cuddel][:user1_public_key]
end

group 'wheel' do
    action  :modify
    members [node[:cuddel][:user1]]
    append  true
end

#########################
#### Directory admin ####
#########################

# These folders are auto-created by Vagrant. If not using Vagrant e.g.
# when using Chef-Solo for provisioning, these folders are explicitly
# created.
dirs = %w{
  log
}

dirs.each do |component|
    the_dir = "/vagrant/#{component}"

    bash 'setup permissions' do
        code <<-EOH
        	# Check if directory exists
            if [ -d #{the_dir} ]
            then
                # Will enter here if the_dir exists,
                echo "#{the_dir} directory exists"
                chmod -R ugo+rwx #{the_dir}
            else
                mkdir -p #{the_dir}
                # chown -R nginx:gigablog-admin #{the_dir}
                chmod -R ugo+rwx #{the_dir}
            fi
        EOH
    end
end


#######################
#### Configure SSH ####
#######################

# Disable root logins and password authentication
bash 'Configure SSH' do
    code <<-EOH
        sed -i -- 's/#PermitRootLogin yes/PermitRootLogin no/g' /etc/ssh/sshd_config
        sed -i -- 's/PasswordAuthentication yes/PasswordAuthentication no/g' /etc/ssh/sshd_config
    EOH
end
