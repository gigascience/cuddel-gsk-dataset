# -*- mode: ruby -*-
# vi: set ft=ruby :

if ENV['CUDDEL_BOX'] == 'aws'
  box = "dummy"
  box_url = "https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box"
else
  box = "centos7"
  box_url = "https://www.dropbox.com/s/0y0uyex2lcxp622/CentOS-7-x86_64-v20170830.box"
end

Vagrant.configure(2) do |config|
  config.vm.box = box
  config.vm.box_url = box_url

  # Cache packages to reduce provisioning time
  if Vagrant.has_plugin?("vagrant-cachier")
    config.cache.scope = :box
  end

  # Forward ports from guest to host
  config.vm.network "forwarded_port", guest: 80, host: 9170
  config.vm.network "forwarded_port", guest: 3306, host: 9171

  # Configure shared folder
  config.vm.synced_folder ".", "/vagrant"

  FileUtils.mkpath("./log")
  FileUtils.chmod_R 0777, ["./log"]

  ####################
  #### VirtualBox ####
  ####################
  config.vm.provider :virtualbox do |vb|
    vb.customize ["setextradata", :id, "VBoxInternal2/SharedFoldersEnableSymlinksCreate//vagrant","1"]
  end

  #############
  #### AWS ####
  #############
  config.vm.provider :aws do |aws, override|
    aws.access_key_id = ENV['AWS_ACCESS_KEY_ID']
    aws.secret_access_key = ENV['AWS_SECRET_ACCESS_KEY']
    aws.keypair_name = ENV['AWS_KEYPAIR_NAME']
    # aws.ami = "ami-1bfa2b78" # selinux disabled
    aws.ami = "ami-b85e86db" # selinux on
    aws.region = ENV['AWS_DEFAULT_REGION']
    aws.instance_type = "t2.micro"
    aws.tags = {
      'Name' => 'GSK analysis',
      'Deployment' => 'test',
    }
    aws.security_groups = ENV['AWS_SECURITY_GROUPS']

    override.ssh.username = "centos"
    override.ssh.private_key_path = ENV['AWS_SSH_PRIVATE_KEY_PATH']
  end

  # Enable provisioning with chef solo
  config.vm.provision :chef_solo do |chef|
    chef.cookbooks_path = [
      "chef/site-cookbooks",
      "chef/chef-cookbooks",
    ]
    chef.environments_path = 'chef/environments'

    ####################################################
    #### Set server environment: development or aws ####
    ####################################################
    chef.environment = "development"

    if ENV['CUDDEL_BOX'] == 'aws'
        chef.add_recipe "aws"
    else
        chef.add_recipe "vagrant"
    end

    # You may also specify custom JSON attributes:
    chef.json = {
      :cuddel_box => ENV['CUDDEL_BOX'],
      :environment => "vagrant"
    }

    # Additional chef settings to put in solo.rb
    chef.custom_config_path = "Vagrantfile.chef"
  end
end
