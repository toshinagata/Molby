function get_download_url() {
  fetch("https://api.github.com/repos/toshinagata/molby/releases")
  .then((response) => {
    console.log(response);
    if (response.ok) {
      return response.json();
    } else {
      console.log("response NG");
    }
  })
  .then((data) => {
    var title0 = data[0].name
    var date0 = new Date(data[0].published_at).toLocaleString();
    document.getElementById("title1").innerHTML = title0 + "<br />" + date0
    for (i = 1; i <= 5; i++) {
      var asset1 = data[0].assets[i - 1];
      var date1 = new Date(asset1.updated_at).toLocaleString();
      var url1 = asset1.browser_download_url;
      var name1 = asset1.name;
      var file1 = "<a href='" + url1 + "'>" + name1 + "</a>";
      var count1 = asset1.download_count;
      document.getElementById("file" + String(i)).innerHTML = file1;
      document.getElementById("count" + String(i)).innerHTML = count1;
    }
    console.log(data);
  })
}
