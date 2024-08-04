#!/usr/bin/ruby
# -*- coding: utf-8 -*-
require "rexml/document"

class REXML::Element
  attr_accessor :up, :down, :prev, :next, :title_en, :title_ja
  def each_element_recursive(&block)
    e = block.call(self)
    return e if e != nil
    i = 1
    while (e = elements[i]) != nil
      ee = e.each_element_recursive(&block)
      if ee != nil && ee != e
        elements[i] = ee
      end
      i += 1
    end
    return nil
  end
end

#  Parse the source document
file = open("| cat src/doc_source*.html", "r")
xml = file.read
file.close
doc = REXML::Document.new(xml)

xmlns = doc.elements["html"].attributes["xmlns"]
head = doc.elements["html/head"]
body = doc.elements["html/body"]

#  Get all 'file' nodes
@file_hash = Hash.new
body.each_element_with_attribute('file') { |e|
  file = e.attributes['file']
  @file_hash[file] = e
  j = 1
  title_en = title_ja = lang = nil
  e.each_element_recursive { |ee|
    if (ln = ee.attributes["lang"]) != nil
      lang = ln
    end
    if ee.name == "h1"
      text = ee.text
      if lang == "ja"
        title_ja = text
      else
        title_en = text
      end
    end
    nil
  }
  e.title_en = title_en
  e.title_ja = title_ja
  while (ee = e.elements[j, "link"]) != nil
    if (rel = ee.attributes["rel"]) != nil
      href = ee.attributes["href"]
      if href == nil || href == ""
        raise "No href attribute at node #{ee.xpath} in #{e.inspect}"
      elsif (el = @file_hash[href]) == nil
        raise "Unknown href '#{href}' at node #{ee.xpath} in #{e.inspect}"
      end
      if rel == "Up"
        e.up = href
        (el.down = (el.down || Array.new)) << file
      elsif rel == "Prev"
        e.prev = href
        el.next = file
      end
    end
    j += 1
  end
}

def replace_special_nodes(e)
  i = 1
  while (ee = elements[i]) != nil
    if ee.name == "link" && (id = ee.attributes["id"]) != nil
      enew = nil
      if id == "\#header"
      elsif id == "\#navigation"
      elsif id == "\#toc"
      end
      if enew != nil
        elements[i] = enew
      end
    end
    i += 1
  end
  if e.name == "link" && (id = e.attributes["id"]) != nil
    if id == "#header"
    end
  end 
end

def make_a_element(href, text)
  a = REXML::Element.new("a")
  a.add_attribute("href", href)
  a.add_text(text)
  a
end

def toc_all(e, level, lang)
  ul = REXML::Element.new("ul")
  e.down.each { |href|
    li = ul.add_element("li")
    ep = @file_hash[href]
    a = make_a_element(href, (lang == "ja" ? ep.title_ja : ep.title_en))
    if level == 0
      li.add_element("h3").add_element(a)
    else
      li.add_element(a)
    end
    if ep.down != nil && ep.down.length > 0
      li.add_element(toc_all(ep, level + 1, lang))
    end
  }
  ul
end

def prev_link(ef)
  return ef.prev if ef.prev
  ef_up = (ef.up && @file_hash[ef.up])
  if ef_up != nil
	ef_up_prev = prev_link(ef_up)
	if ef_up_prev && (down = @file_hash[ef_up_prev].down) && down.length > 0
	  return down[-1]
	end
  end
  return nil
end

def next_link(ef)
  return ef.next if ef.next
  ef_up = (ef.up && @file_hash[ef.up])
  if ef_up != nil
	ef_up_next = next_link(ef_up)
	if ef_up_next && (down = @file_hash[ef_up_next].down) && down.length > 0
	  return down[0]
	end
  end
  return nil
end
	  
def special_node(e, ef, lang)
  return nil if e.name != "link" || (id = e.attributes["id"]) == nil
  if id == "\#header"
    en = REXML::Element.new("div")
    en.add_attribute("class", "topic_path")
    ep = ef
    c = [(lang == "ja" ? ef.title_ja : ef.title_en)]
    while ep.up != nil
      c.unshift(REXML::Text.new(" &gt; ", false, nil, true))
      href = ep.up
      ep = @file_hash[href]
      c.unshift(make_a_element(href, (lang == "ja" ? ep.title_ja : ep.title_en)))
    end
    if (href = prev_link(ef)) != nil
      c.push(REXML::Text.new(" &nbsp;&nbsp; ", false, nil, true))
      c.push(make_a_element(href, (lang == "ja" ? "[前]" : "[Prev]")))
    end
    if (href = next_link(ef)) != nil
      c.push(REXML::Text.new(" &nbsp;&nbsp; ", false, nil, true))
      c.push(make_a_element(href, (lang == "ja" ? "[次]" : "[Next]")))
    end
    otherlang = (lang == "ja" ? "en" : "ja")
    n = REXML::Element.new("span")
    n.add_attribute("class", "float_right")
    n.add(make_a_element('../' + otherlang + '/' + ef.attributes['file'], (otherlang == "ja" ? "[Japanese]" : "[English]")))
    c.push(n)
    c.each { |n|
      if n.is_a?(String)
        en.add_text(n)
      else
        en.add(n)
      end
    }
  elsif id == "\#navigation"
    en = REXML::Element.new("div")
    en.add_attribute("class", "navigation")
    if (href = ef.up) != nil
      if href != "index.html"
        en.add(make_a_element("index.html", (lang == "ja" ? "[トップ]" : "[Top]")))
        en.add_text(" ")
      end
      ep = @file_hash[href]
      en.add(make_a_element(href, (lang == "ja" ? "[上: #{ep.title_ja}]" : "[Up: #{ep.title_en}]")))
      en.add_text(" ")
    else
      en.add_text(lang == "ja" ? "[トップ] " : "[Top] ")
    end
    if (href = prev_link(ef)) != nil
      ep = @file_hash[href]
      en.add(make_a_element(href, (lang == "ja" ? "[前: #{ep.title_ja}]" : "[Prev: #{ep.title_en}]")))
      en.add_text(" ")
    end
    if (href = next_link(ef)) != nil
      ep = @file_hash[href]
      en.add(make_a_element(href, (lang == "ja" ? "[次: #{ep.title_ja}]" : "[Next: #{ep.title_en}]")))
      en.add_text(" ")
    end
    otherlang = (lang == "ja" ? "en" : "ja")
    n = REXML::Element.new("span")
    n.add_attribute("class", "float_right")
    n.add(make_a_element('../' + otherlang + '/' + ef.attributes['file'], (otherlang == "ja" ? "[Japanese]" : "[English]")))
    en.add(n)
  elsif id == "\#toc"
    en = REXML::Element.new("div")
    en.add_attribute("class", "contents")
	if ef.down
      ul = en.add_element("ul")
      ef.down.each { |href|
        li = ul.add_element("li")
        ep = @file_hash[href]
        li.add_element(make_a_element(href, (lang == "ja" ? ep.title_ja : ep.title_en)))
      }
	end
  elsif id == "\#toc_all"
    en = REXML::Element.new("div")
    en.add_attribute("class", "contents")
    en.add_element(toc_all(ef, 0, lang))
  else
    en = nil
  end
  en
end

base_dir = "../docs/MolbyDoc"
system("mkdir -p #{base_dir}; rm -rf #{base_dir}")

#  Output to files (en and jp)
preamble = doc.xml_decl.to_s + "\n" + doc.doctype.to_s + "\n"
body.each_element_with_attribute('file') { |ef|
  file = ef.attributes['file']
  for lang in ["ja", "en"]
    system("mkdir -p #{base_dir}/#{lang}")
    ndoc = REXML::Document.new(preamble)
    #  Root element (html)
    html = REXML::Element.new("html")
	html.attributes["lang"] = lang
	html.attributes["xml:lang"] = lang
	if xmlns
	  html.attributes["xmlns"] = xmlns
	end
    ndoc.add(html)
    #  head section
    nhead = html.add_element("head")
    head.each_element { |e|
      e = e.deep_clone
      if e.name == "title"
        e = REXML::Element.new("title")
        e.add_text(lang == "ja" ? ef.title_ja : ef.title_en)
      elsif e.name == "meta"
        if e.attributes["http-equiv"] == "Content-Type"
          e.attributes["lang"] = lang
		  e.attributes["xml:lang"] = lang
        end
      end
      nhead.add_element(e)
    }
    nhead.add_element("link", "rel"=>"Start", "href"=>"index.html")
    if ef.up != nil && ef.up != ""
      nhead.add_element("link", "rel"=>"Index", "href"=>ef.up)
    end
    if ef.prev != nil && ef.prev != ""
      nhead.add_element("link", "rel"=>"Prev", "href"=>ef.prev)
    end
    if ef.next != nil && ef.next != ""
      nhead.add_element("link", "rel"=>"Next", "href"=>ef.next)
    end
    #  body section
    bd = html.add_element("body")
    ef.each_element { |e|
      next if (la = e.attributes["lang"]) != nil && la != lang
      next if e.name == "link" && e.attributes["id"] == nil
	  if la
	    e.attributes["xml:lang"] = la
	  end
      bd.add_element(e.deep_clone)
    }
    #  Convert special nodes
    bd.each_element_recursive { |e|
      if e.name == "link" && e.attributes["id"] != nil
        val = special_node(e, ef, lang)
      else
	    if e.name == "div" && e.attributes["lang"] != nil
		  e.attributes["xml:lang"] = e.attributes["lang"]
		end
        nil
      end
    }
    open("#{base_dir}/#{lang}/_#{file}", "w") { |fp|
      ndoc.write(fp)
      fp.write("\n")
    }
	system("cd #{base_dir}/#{lang}; sed 's!/>! />!g' _#{file} >#{file}; rm _#{file}")
  end
}

system("cp -r -p etc #{base_dir}; rm -rf #{base_dir}/etc/CVS #{base_dir}/etc/.svn #{base_dir}/etc/.DS_Store") 
system("cp -r -p src/molby_rb #{base_dir}/en; rm -rf #{base_dir}/en/molby_rb/CVS #{base_dir}/en/molby_rb/.svn")
system("cp -r -p #{base_dir}/en/molby_rb #{base_dir}/ja")
print "Documents were successfully created in #{base_dir}/{en,ja}.\n"
exit 0
