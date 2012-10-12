#
#  dialog.rb
#
#  Created by Toshi Nagata on 2008/06/27.
#  Copyright 2008 Toshi Nagata. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

class Dialog

  def self.run(*args, &block)
    obj = Dialog.new(*args)
    obj.instance_eval(&block)
    obj.run
  end

  alias initialize_orig initialize

  def initialize(*args, &block)
    initialize_orig(*args)
	instance_eval(&block) if block
  end

  def value(tag)
    attr(tag, :value)
  end
  
  def set_value(tag, value)
    set_attr(tag, :value=>value)
	value
  end

  def self.filter_kit(title, description, &block)
    Dialog.new(title, nil, nil) {
      def write(s)  #  Override standard output
        item_with_tag("text").append_string(s)
      end
      button_action = proc {
        names = Dialog.open_panel("Select file(s) to process", nil, nil, false, true)
        if names
          stdout_save = $stdout
          $stdout = self
          block.call(names)
          $stdout = stdout_save
        end
      }
      layout(1,
        item(:text, :title=>description),
        item(:button, :title=>"Select Files...", :action=>button_action),
        item(:textview, :width=>320, :height=>200, :editable=>false, :tag=>"text", :font=>[:fixed, 10]),
        item(:button, :title=>"Exit", :action=>proc { hide }, :align=>:center))
      show
    }
    nil
  end
  
end
